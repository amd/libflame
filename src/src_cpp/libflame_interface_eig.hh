/******************************************************************************
 * Copyright (C) 2021-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface_eig.hh
 *  libflame_interface.hh defines all the eigenvalue routines for Libflame CPP templated public
 *  interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_EIG_HH
#define LIBFLAME_INTERFACE_EIG_HH

#include "libflame.hh"

namespace libflame
{
     /** @defgroup Eig Eigenvalue Routines
     * @ingroup LAPACK
     * @{
     */

    /** @defgroup NYS Non-symmetric eigenvalues
     * @ingroup Eig
     * @{
     */

    /** @defgroup geev geev
     * @ingroup NYS
     * @{
     */
    /*! @brief GEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for
GE matrices
*
* @details
* \b Purpose:
* \verbtatim
 GEEV computes for an N-by-N real nonsymmetric matrix A, the
 eigenvalues and, optionally, the left and/or right eigenvectors.

 The right eigenvector v(j) of A satisfies
                  A * v(j) = lambda(j) * v(j)
 where lambda(j) is its eigenvalue.
 The left eigenvector u(j) of A satisfies
               u(j)**H * A = lambda(j) * u(j)**H
 where u(j)**H denotes the conjugate-transpose of u(j).

 The computed eigenvectors are normalized to have Euclidean norm
 equal to 1 and largest component real.
  \endverbatim

 * @param[in] JOBVL
          JOBVL is CHARACTER*1 \n
          = 'N': left eigenvectors of A are not computed; \n
          = 'V': left eigenvectors of A are computed. \n
 * @param[in] JOBVR
          JOBVR is CHARACTER*1 \n
          = 'N': right eigenvectors of A are not computed; \n
          = 'V': right eigenvectors of A are computed. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
          On exit, A has been overwritten. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] WR
          WR is REAL array, dimension (N) \n
 * @param[out] WI
          WI is REAL array, dimension (N) \n
          WR and WI contain the real and imaginary parts,
          respectively, of the computed eigenvalues.  Complex
          conjugate pairs of eigenvalues appear consecutively
          with the eigenvalue having the positive imaginary part
          first. \n
 * @param[out] VL
          VL is REAL array, dimension (LDVL,N) \n
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order
          as their eigenvalues. \n
          If JOBVL = 'N', VL is not referenced. \n
          If the j-th eigenvalue is real, then u(j) = VL(:,j),
          the j-th column of VL. \n
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair - computes row and column scaling to reduce condition number of matrixhen
u(j@ = VL(:,j) + i*VL(:,j+1) and u(j+1) = VL(:,j) \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the array VL.  LDVL >= 1; if
          JOBVL = 'V', LDVL >= N. \n
 * @param[out] VR
          VR is REAL array, dimension (LDVR,N) \n
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order
          as their eigenvalues. \n
          If JOBVR = 'N', VR is not referenced. \n
          If the j-th eigenvalue is real, then v(j) = VR(:,j),
          the j-th column of VR. \n
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
          v(j+1) = VR(:,j) - i*VR(:,j+1). \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR.  LDVR >= 1; if
          JOBVR = 'V', LDVR >= N. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,2*N).
          For good performance, LWORK must generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, the QR algorithm failed to compute all the
                eigenvalues, and no eigenvectors have been computed;
                elements i+1:N of W contain eigenvalues which have
                converged. \n

 *  * */
    template <typename T>
    void geev(char *jobvl, char *jobvr, integer *n, T *a, integer *lda, T *wr, T *wi, T *vl,
              integer *ldvl, T *vr, integer *ldvr, T *work, integer *lwork, integer *info)
    {
        geev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
    }
    template <typename T, typename Ta>
    void geev(char *jobvl, char *jobvr, integer *n, T *a, integer *lda, T *w, T *vl, integer *ldvl,
              T *vr, integer *ldvr, T *work, integer *lwork, Ta *rwork, integer *info)
    {
        geev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
    }
    /** @}*/ // end of geev

    /** @defgroup geevx geevx
     * @ingroup NYS
     * @{
     */
    /*! @brief GEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors
for GE matrices(enabling conditions)
* @details
* \b Purpose:
* \verbatim
    GEEVX computes for an N-by-N real nonsymmetric matrix A, the
    eigenvalues and, optionally, the left and/or right eigenvectors.

    Optionally also, it computes a balancing transformation to improve
    the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
    SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
    (RCONDE), and reciprocal condition numbers for the right
    eigenvectors (RCONDV).

    The right eigenvector v(j) of A satisfies
                    A * v(j) = lambda(j) * v(j)
    where lambda(j) is its eigenvalue.
    The left eigenvector u(j) of A satisfies
                 u(j)**H * A = lambda(j) * u(j)**H
    where u(j)**H denotes the conjugate-transpose of u(j).

    The computed eigenvectors are normalized to have Euclidean norm
    equal to 1 and largest component real.

    Balancing a matrix means permuting the rows and columns to make it
    more nearly upper triangular, and applying a diagonal similarity
    transformation D * A * D**(-1), where D is a diagonal matrix, to
    make its rows and columns closer in norm and the condition numbers
    of its eigenvalues and eigenvectors smaller.  The computed
    reciprocal condition numbers correspond to the balanced matrix.
    Permuting rows and columns will not change the condition numbers
    (in exact arithmetic) but diagonal scaling will.  For further
    explanation of balancing, see section 4.10.2 of the LAPACK
    Users' Guide.
 \endverbatim

 * @param[in] BALANC
          BALANC is CHARACTER*1 \n
          Indicates how the input matrix should be diagonally scaled
          and/or permuted to improve the conditioning of its
          eigenvalues. \n
          = 'N': Do not diagonally scale or permute; \n
          = 'P': Perform permutations to make the matrix more nearly
                 upper triangular. Do not diagonally scale; \n
          = 'S': Diagonally scale the matrix, i.e. replace A by
                 D*A*D**(-1), where D is a diagonal matrix chosen
                 to make the rows and columns of A more equal in
                 norm. Do not permute; \n
          = 'B': Both diagonally scale and permute A. \n
 \n
          Computed reciprocal condition numbers will be for the matrix
          after balancing and/or permuting. Permuting does not change
          condition numbers (in exact arithmetic), but balancing does. \n
 * @param[in] JOBVL
          JOBVL is CHARACTER*1 \n
          = 'N': left eigenvectors of A are not computed; \n
          = 'V': left eigenvectors of A are computed. \n
          If SENSE = 'E' or 'B', JOBVL must = 'V'. \n
 * @param[in] JOBVR
          JOBVR is CHARACTER*1 \n
          = 'N': right eigenvectors of A are not computed; \n
          = 'V': right eigenvectors of A are computed. \n
          If SENSE = 'E' or 'B', JOBVR must = 'V'. \n
 * @param[in] SENSE
          SENSE is CHARACTER*1 \n
          Determines which reciprocal condition numbers are computed. \n
          = 'N': None are computed; \n
          = 'E': Computed for eigenvalues only; \n
          = 'V': Computed for right eigenvectors only; \n
          = 'B': Computed for eigenvalues and right eigenvectors. \n
 \n
          If SENSE = 'E' or 'B', both left and right eigenvectors
          must also be computed (JOBVL = 'V' and JOBVR = 'V'). \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
          On exit, A has been overwritten.  If JOBVL = 'V' or
          JOBVR = 'V', A contains the real Schur form of the balanced
          version of the input matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] WR
          WR is REAL array, dimension (N) \n
 * @param[out] WI
          WI is REAL array, dimension (N) \n
          WR and WI contain the real and imaginary parts,
          respectively, of the computed eigenvalues.  Complex
          conjugate pairs of eigenvalues will appear consecutively
          with the eigenvalue having the positive imaginary part
          first. \n
 * @param[out] VL
          VL is REAL array, dimension (LDVL,N) \n
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order
          as their eigenvalues. \n
          If JOBVL = 'N', VL is not referenced. \n
          If the j-th eigenvalue is real, then u(j) = VL(:,j),
          the j-th column of VL. \n
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
          u(j+1) = VL(:,j) - i*VL(:,j+1). \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the array VL.  LDVL >= 1; if
          JOBVL = 'V', LDVL >= N. \n
 * @param[out] VR
          VR is REAL array, dimension (LDVR,N) \n
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order
          as their eigenvalues. \n
          If JOBVR = 'N', VR is not referenced. \n
          If the j-th eigenvalue is real, then v(j) = VR(:,j),
          the j-th column of VR. \n
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
          v(j+1) = VR(:,j) - i*VR(:,j+1). \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR.  LDVR >= 1, and if
          JOBVR = 'V', LDVR >= N. \n
 * @param[out] ILO
          ILO is INTEGER \n
 * @param[out] IHI
          IHI is INTEGER \n
          ILO and IHI are integer values determined when A was
          balanced.  The balanced A(i,j) = 0 if I > J and
          J = 1,...,ILO-1 or I = IHI+1,...,N. \n
 * @param[out] SCALE
          SCALE is REAL array, dimension (N) \n
          Details of the permutations and scaling factors applied
          when balancin - computes row and column scaling to reduce condition number of matrix.  If
P@j) is the index of the row and column interchanged with row and column j, and D(j) is the
scaling factor applied to row and column j, then SCALE(J) = P(J),    for J = 1,...,ILO-1 = D(J), for
J = ILO,...,IHI = P(J)     for J = IHI+1,...,N. The order in which the interchanges are made is N to
IHI+1, then 1 to ILO-1. \n
 * @param[out] ABNRM
          ABNRM is REAL \n
          The one-norm of the balanced matrix (the maximum
          of the sum of absolute values of elements of any column). \n
 * @param[out] RCONDE
          RCONDE is REAL array, dimension (N) \n
          RCONDE(j) is the reciprocal condition number of the j-th
          eigenvalue. \n
 * @param[out] RCONDV
          RCONDV is REAL array, dimension (N) \n
          RCONDV(j) is the reciprocal condition number of the j-th
          right eigenvector. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  If SENSE = 'N' or 'E',
          LWORK >= fla_max(1,2*N), and if SENSE = 'V' or 'B',
          LWORK >= N*N+2*N. \n
          For good performance, LWORK must generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, the QR algorithm failed to compute all the
                eigenvalues, and no eigenvectors or condition numbers
                have been computed; elements 1:ILO-1 and i+1:N of W
                contain eigenvalues which have converged. \n

 *  * */
    template <typename T>
    void geevx(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, T *a, integer *lda,
               T *wr, T *wi, T *vl, integer *ldvl, T *vr, integer *ldvr, integer *ilo, integer *ihi,
               T *scale, T *abnrm, T *rconde, T *rcondv, T *work, integer *lwork, integer *iwork,
               integer *info)
    {
        geevx(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale,
              abnrm, rconde, rcondv, work, lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void geevx(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, T *a, integer *lda,
               T *w, T *vl, integer *ldvl, T *vr, integer *ldvr, integer *ilo, integer *ihi,
               Ta *scale, Ta *abnrm, Ta *rconde, Ta *rcondv, T *work, integer *lwork, Ta *rwork,
               integer *info)
    {
        geevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm,
              rconde, rcondv, work, lwork, rwork, info);
    }
    /** @}*/ // end of geevx

    /** @defgroup gees gees
     * @ingroup NYS
     * @{
     */
    /*! @brief GEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
 vectors for GE matrices</b>
 *  @details
 *  \b Purpose:
 *  \verbatim
    GEES computes for an N-by-N real nonsymmetric matrix A, the
    eigenvalues, the real Schur form T, and, optionally, the matrix of
    Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).

    Optionally, it also orders the eigenvalues on the diagonal of the
    real Schur form so that selected eigenvalues are at the top left.
    The leading columns of Z then form an orthonormal basis for the
    invariant subspace corresponding to the selected eigenvalues.

    A matrix is in real Schur form if it is upper quasi-triangular with
    1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
    form
           [  a  b  ]
           [  c  a  ]

    where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
    \endverbatim

 * @param[in] JOBVS
          JOBVS is CHARACTER*1 \n
          = 'N': Schur vectors are not computed; \n
          = 'V': Schur vectors are computed. \n
 * @param[in] SORT
          SORT is CHARACTER*1 \n
          Specifies whether or not to order the eigenvalues on the
          diagonal of the Schur form. \n
          = 'N': Eigenvalues are not ordered; \n
          = 'S': Eigenvalues are ordered (see SELECT). \n
 * @param[in] SELECT
          SELECT is a LOGICAL FUNCTION of two REAL arguments \n
          SELECT must be declared EXTERNAL in the calling subroutine.
          If SORT = 'S', SELECT is used to select eigenvalues to sort
          to the top left of the Schur form.
          If SORT = 'N', SELECT is not referenced. \n
          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
          conjugate pair of eigenvalues is selected, then both complex
          eigenvalues are selected. \n
          Note that a selected complex eigenvalue may no longer
          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
          ordering may change the value of complex eigenvalues
          (especially if the eigenvalue is ill-conditioned); in this
          case INFO is set to N+2 (see INFO below). \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
          On exit, A has been overwritten by its real Schur form T. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] SDIM
          SDIM is INTEGER \n
          If SORT = 'N', SDIM = 0. \n
          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
                         for which SELECT is true. (Complex conjugate
                         pairs for which SELECT is true for either
                         eigenvalue count as  - computes row and column scaling to reduce condition
 number of matrix \n
 * @param[out] WR
          WR is REAL array, dimension (N) \n
 * @param[out] WI
          WI is REAL array, dimension (N) \n
          WR and WI contain the real and imaginary parts,
          respectively, of the computed eigenvalues in the same order
          that they appear on the diagonal of the output Schur form T.
          Complex conjugate pairs of eigenvalues will appear
          consecutively with the eigenvalue having the positive
          imaginary part first. \n
 * @param[out] VS
          VS is REAL array, dimension (LDVS,N) \n
          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
          vectors. \n
          If JOBVS = 'N', VS is not referenced. \n
 * @param[in] LDVS
          LDVS is INTEGER \n
          The leading dimension of the array VS.  LDVS >= 1; if
          JOBVS = 'V', LDVS >= N. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,2*N).
          For good performance, LWORK must generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (N) \n
 * @param[out]	BWORK
          BWORK is LOGICAL array, dimension (N) \n
          Not referenced if SORT = 'N'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value. \n
          > 0: if INFO = i, and i is \n
               <= N:  the QR algorithm failed to compute all the
                      eigenvalues; elements 1:ILO-1 and i+1:N of W
                      contain those eigenvalues which have converged;
                      if JOBVS = 'V', VS contains the matrix which
                      reduces A to its partially converged Schur form. \n
               = N+1: the eigenvalues could not be reordered because
                      some eigenvalues were too close to separate (the
                      problem is very ill-conditioned); \n
               = N+2: after reordering, roundoff changed values of
                      some complex eigenvalues so that leading
                      eigenvalues in the Schur form no longer satisfy
                      SELECT = .TRUE..  This could also be caused by
                      underflow due to scaling. \n

 *  * */
    template <typename T>
    void gees(char *jobvs, char *sort, void *select, integer *n, T *a, integer *lda, integer *sdim,
              T *wr, T *wi, T *vs, integer *ldvs, T *work, integer *lwork, logical *bwork,
              integer *info)
    {
        gees(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info);
    }
    template <typename T, typename Ta>
    void gees(char *jobvs, char *sort, void *select, integer *n, T *a, integer *lda, integer *sdim,
              T *w, T *vs, integer *ldvs, T *work, integer *lwork, Ta *rwork, logical *bwork,
              integer *info)
    {
        gees(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info);
    }
    /** @}*/ // end of gees

    /** @defgroup geesx geesx
     * @ingroup NYS
     * @{
     */
    /*! @brief GEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
 vectors for GE matrices
 *
 *  @details
 *  \b Purpose:
 *  \verbatim
    GEESX computes for an N-by-N real nonsymmetric matrix A, the
    eigenvalues, the real Schur form T, and, optionally, the matrix of
    Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).

    Optionally, it also orders the eigenvalues on the diagonal of the
    real Schur form so that selected eigenvalues are at the top left;
    computes a reciprocal condition number for the average of the
    selected eigenvalues (RCONDE); and computes a reciprocal condition
    number for the right invariant subspace corresponding to the
    selected eigenvalues (RCONDV).  The leading columns of Z form an
    orthonormal basis for this invariant subspace.

    For further explanation of the reciprocal condition numbers RCONDE
    and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where
    these quantities are called s and sep respectively).

    A real matrix is in real Schur form if it is upper quasi-triangular
    with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in
    the form
             [  a  b  ]
             [  c  a  ]

    where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
    \endverbatim

 * @param[in] JOBVS
          JOBVS is CHARACTER*1 \n
          = 'N': Schur vectors are not computed; \n
          = 'V': Schur vectors are computed. \n
 * @param[in] SORT
          SORT is CHARACTER*1 \n
          Specifies whether or not to order the eigenvalues on the
          diagonal of the Schur form. \n
          = 'N': Eigenvalues are not ordered; \n
          = 'S': Eigenvalues are ordered (see SELECT). \n
 * @param[in] SELECT
          SELECT is a LOGICAL FUNCTION of two REAL arguments
          SELECT must be declared EXTERNAL in the calling subroutine. \n
          If SORT = 'S', SELECT is used to select eigenvalues to sort
          to the top left of the Schur form.
          If SORT = 'N', SELECT is not referenced. \n
          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
          SELECT(WR(j),WI(j)) is true; i.e., if either one of a
          complex conjugate pair of eigenvalues is selected, then both
          are.  Note that a selected complex eigenvalue may no longer
          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
          ordering may change the value of complex eigenvalues
          (especially if the eigenvalue is ill-conditioned); in this
          case INFO may be set to N+3 (see INFO below). \n
 * @param[in] SENSE
          SENSE is CHARACTER*1 \n
          Determines which reciprocal condition numbers are computed. \n
          = 'N': None are computed; \n
          = 'E': Computed for average of selected eigenvalues only; \n
          = 'V': Computed for selected right invariant subspace only; \n
          = 'B': Computed for both. \n
          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the N-by-N matrix A. \n
          On exit, A is overwritten by its real Schur form T. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] SDIM
          SDIM is INTEGER \n
          If SORT = 'N', SDIM = 0. \n
          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
                         for which SELECT is true. (Complex conjugate
                         pairs for which SELECT is true for either
                         eigenvalue count as 2.) \n
 * @param[out] WR
          WR is REAL array, dimension (N) \n
 * @param[out] WI
          WI is REAL array, dimension (N) \n
          WR and WI contain the real and imaginary parts, respectively,
          of the computed eigenvalues, in the same order that they
          appear on the diagonal of the output Schur form T.  Complex
          conjugate pairs of eigenvalues appear consecutively with the
          eigenvalue having the positive imaginary part fir -
          computes row and column scaling to reduce condition number of matrix \n
 * @param[out] VS
          VS is REAL array, dimension (LDVS,N) \n
          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
          vectors. \n
          If JOBVS = 'N', VS is not referenced. \n
 * @param[in] LDVS
          LDVS is INTEGER \n
          The leading dimension of the array VS.  LDVS >= 1, and if
          JOBVS = 'V', LDVS >= N. \n
 * @param[out] RCONDE
          RCONDE is REAL \n
          If SENSE = 'E' or 'B', RCONDE contains the reciprocal
          condition number for the average of the selected eigenvalues.
          Not referenced if SENSE = 'N' or 'V'. \n
 * @param[out] RCONDV
          RCONDV is REAL \n
          If SENSE = 'V' or 'B', RCONDV contains the reciprocal
          condition number for the selected right invariant subspace.
          Not referenced if SENSE = 'N' or 'E'. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,2*N). \n
          Also, if SENSE = 'E' or 'V' or 'B', LWORK >= 2*SDIM*(N-SDIM),
          where SDIM is the number of selected eigenvalues computed by
          this routine.  Note that 2*SDIM*(N-SDIM) <= N*N/2. Note also
          that an error is only returned if LWORK < fla_max(1,2*N), but if
          SENSE = 'E' or 'V' or 'B' this may not be large enough.
          For good performance, LWORK must generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates upper bound on the optimal size of the
          array WORK, returns this value as the first entry of the WORK
          array, and no error message related to LWORK is issued by
          XERBLA. \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (N) \n
 * @param[out]	BWORK
          BWORK is LOGICAL array, dimension (N) \n
          Not referenced if SORT = 'N'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value. \n
          > 0: if INFO = i, and i is \n
             <= N: the QR algorithm failed to compute all the
                   eigenvalues; elements 1:ILO-1 and i+1:N of W
                   contain those eigenvalues which have converged; if
                   JOBVS = 'V', VS contains the transformation which
                   reduces A to its partially converged Schur form. \n
             = N+1: the eigenvalues could not be reordered because some
                   eigenvalues were too close to separate (the problem
                   is very ill-conditioned); \n
             = N+2: after reordering, roundoff changed values of some
                   complex eigenvalues so that leading eigenvalues in
                   the Schur form no longer satisfy SELECT=.TRUE.  This
                   could also be caused by underflow due to scaling. \n

 *  * */
    template <typename T>
    void geesx(char *jobvs, char *sort, void *select, char *sense, integer *n, T *a, integer *lda,
               integer *sdim, T *wr, T *wi, T *vs, integer *ldvs, T *rconde, T *rcondv, T *work,
               integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info)
    {
        geesx(jobvs, sort, select, sense, n, a, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work,
              lwork, iwork, liwork, bwork, info);
    }
    template <typename T, typename Ta>
    void geesx(char *jobvs, char *sort, void *select, char *sense, integer *n, T *a, integer *lda,
               integer *sdim, T *w, T *vs, integer *ldvs, Ta *rconde, Ta *rcondv, T *work,
               integer *lwork, Ta *rwork, logical *bwork, integer *info)
    {
        geesx(jobvs, sort, select, sense, n, a, lda, sdim, w, vs, ldvs, rconde, rcondv, work, lwork,
              rwork, bwork, info);
    }
    /** @}*/ // end of geesx

    /** @defgroup ggev3 ggev3
     * @ingroup NYS
     * @{
     */
    /*! @brief GGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors
for GE matrices (blocked algorithm)
 * @details
 * \b Purpose:
    \verbatim
     GGEV3 computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     the generalized eigenvalues, and optionally, the left and/or right
     generalized eigenvectors.
     A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     singular. It is usually represented as the pair (alpha,beta), as
     there is a reasonable interpretation for beta=0, and even for both
     being zero.
     The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     of (A,B) satisfies
                      A * v(j) = lambda(j) * B * v(j).
     The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     of (A,B) satisfies
                      u(j)**H * A  = lambda(j) * u(j)**H * B .
     where u(j)**H is the conjugate-transpose of u(j).
    \endverbatim

 * @param[in] JOBVL
          JOBVL is CHARACTER*1 \n
          = 'N':  do not compute the left generalized eigenvectors; \n
          = 'V':  compute the left generalized eigenvectors. \n
 * @param[in] JOBVR
          JOBVR is CHARACTER*1 \n
          = 'N':  do not compute the right generalized eigenvectors; \n
          = 'V':  compute the right generalized eigenvectors. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, VL, and VR.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the matrix A in the pair (A,B).
          On exit, A has been overwritten. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the matrix B in the pair (A,B).
          On exit, B has been overwritten. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B.  LDB >= fla_max(1,N). \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
          the j-th eigenvalue is real; if positive, then the j-th and
          (j+1)-st eigenvalues are a complex conjugate pair, with
          ALPHAI(j+1) negative. \n
 \n
          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
          may easily over- or underflow, and BETA(j) may even be zero.
          Thus, the user should avoid naively computing the ratio
          alpha/beta.  However, ALPHAR and ALPHAI will be always less
          than and usually comparable with norm(A) in magnitude, and
          BETA always less than and usually comparable with norm(B). \n
 * @param[out] VL
          VL is REAL array, dimension (LDVL,N) \n
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order as
          their eigenvalues. If the j-th eigenvalue is real, then
          u(j) = VL(:,j), the j-th column of VL. If the j-th and
          (j+1)-th eigenvalues form a complex conjugate pair, then
          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
          Each eigenvector is scaled so the largest component has
          abs(real part)+abs(imag. part)=1. \n
          Not referenced if JOBVL = 'N'. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the matrix VL. LDVL >= 1, and
          if JOBVL = 'V', LDVL >= N. \n
 * @param[out] VR
          VR is REAL array, dimension (LDVR,N) \n
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order as
          their eigenvalues. If the j-th eigenvalue is real, then
          v(j) = VR(:,j), the j-th column of VR. If the j-th and
          (j+1)-th eigenvalues form a complex conjugate pair, then
          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
          Each eigenvector is scaled so the largest component has
          abs(real part)+abs(imag. part)=1. \n
          Not referenced if JOBVR = 'N'. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the matrix VR. LDVR >= 1, and
          if JOBVR = 'V', LDVR >= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1,...,N: \n
                The QZ iteration failed.  No eigenvectors have been
                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
                should be correct for j=INFO+1,...,N. \n
          > N:  =N+1: other than QZ iteration failed in SHGEQZ. \n
                =N+2: error return from STGEVC. \n

*  * */
    template <typename T>
    void ggev3(char *jobvl, char *jobvr, integer *n, T *a, integer *lda, T *b, integer *ldb,
               T *alphar, T *alphai, T *beta, T *vl, integer *ldvl, T *vr, integer *ldvr, T *work,
               integer *lwork, integer *info)
    {
        ggev3(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work,
              lwork, info);
    }
    template <typename T, typename Ta>
    void ggev3(char *jobvl, char *jobvr, integer *n, T *a, integer *lda, T *b, integer *ldb,
               T *alpha, T *beta, T *vl, integer *ldvl, T *vr, integer *ldvr, T *work,
               integer *lwork, Ta *rwork, integer *info)
    {
        ggev3(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork,
              info);
    }
    /** @}*/ // end of gg

    /** @defgroup ggev ggev
     * @ingroup NYS
     * @{
     */
    /*! @brief GGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 GE matrices
 * @details
 * \b Purpose:
    \verbatim
     GGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     the generalized eigenvalues, and optionally, the left and/or right
     generalized eigenvectors.
     A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     singular. It is usually represented as the pair (alpha,beta), as
     there is a reasonable interpretation for beta=0, and even for both
     being zero.
     The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     of (A,B) satisfies
                      A * v(j) = lambda(j) * B * v(j).
     The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     of (A,B) satisfies
                      u(j)**H * A  = lambda(j) * u(j)**H * B .
     where u(j)**H is the conjugate-transpose of u(j).
    \endverbatim

 * @param[in] JOBVL
          JOBVL is CHARACTER*1 \n
          = 'N':  do not compute the left generalized eigenvectors; \n
          = 'V':  compute the left generalized eigenvectors. \n
 * @param[in] JOBVR
          JOBVR is CHARACTER*1 \n
          = 'N':  do not compute the right generalized eigenvectors; \n
          = 'V':  compute the right generalized eigenvectors. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, VL, and VR.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the matrix A in the pair (A,B).
          On exit, A has been overwritten. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the matrix B in the pair (A,B).
          On exit, B has been overwritten. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B.  LDB >= fla_max(1,N). \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
          the j-th eigenvalue is real; if positive, then the j-th and
          (j+1)-st eigenvalues are a complex conjugate pair, with
          ALPHAI(j+1) negative. \n
 \n
          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
          may easily over- or underflow, and BETA(j) may even be zero.
          Thus, the user should avoid naively computing the ratio
          alpha/beta.  However, ALPHAR and ALPHAI will be always less
          than and usually comparable with norm(A) in magnitude, and
          BETA always less than and usually comparable with norm(B). \n
 * @param[out] VL
          VL is REAL array, dimension (LDVL,N) \n
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order as
          their eigenvalues. If the j-th eigenvalue is real, then
          u(j) = VL(:,j), the j-th column of VL. If the j-th and
          (j+1)-th eigenvalues form a complex conjugate pair, then
          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
          Each eigenvector is scaled so the largest component has
          abs(real part)+abs(imag. part)=1. \n
          Not referenced if JOBVL = 'N'. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the matrix VL. LDVL >= 1, and
          if JOBVL = 'V', LDVL >= N. \n
 * @param[out] VR
          VR is REAL array, dimension (LDVR,N) \n
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order as
          their eigenvalues. If the j-th eigenvalue is real, then
          v(j) = VR(:,j), the j-th column of VR. If the j-th and
          (j+1)-th eigenvalues form a complex conjugate pair, then
          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
          Each eigenvector is scaled so the largest component has
          abs(real part)+abs(imag. part)=1. \n
          Not referenced if JOBVR = 'N'. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the matrix VR. LDVR >= 1, and
          if JOBVR = 'V', LDVR >= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,8*N).
          For good performance, LWORK must generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1,...,N: \n
                The QZ iteration failed.  No eigenvectors have been
                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
                should be correct for j=INFO+1,...,N. \n
          > N:  =N+1: other than QZ iteration failed in SHGEQZ.
                =N+2: error return from STGEVC. \n

 *  * */
    template <typename T>
    void ggev(char *jobvl, char *jobvr, integer *n, T *a, integer *lda, T *b, integer *ldb,
              T *alphar, T *alphai, T *beta, T *vl, integer *ldvl, T *vr, integer *ldvr, T *work,
              integer *lwork, integer *info)
    {
        ggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork,
             info);
    }
    template <typename T, typename Ta>
    void ggev(char *jobvl, char *jobvr, integer *n, T *a, integer *lda, T *b, integer *ldb,
              T *alpha, T *beta, T *vl, integer *ldvl, T *vr, integer *ldvr, T *work,
              integer *lwork, Ta *rwork, integer *info)
    {
        ggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork,
             info);
    }
    /** @}*/ // end of ggev

    /** @defgroup ggevx ggevx
     * @ingroup NYS
     * @{
     */
    /*! @brief GGEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors
 for GE matrices
 * @details
 * \b Purpose:
    \verbatim
     SGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     the generalized eigenvalues, and optionally, the left and/or right
     generalized eigenvectors.
     Optionally also, it computes a balancing transformation to improve
     the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     the eigenvalues (RCONDE), and reciprocal condition numbers for the
     right eigenvectors (RCONDV).
     A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     singular. It is usually represented as the pair (alpha,beta), as
     there is a reasonable interpretation for beta=0, and even for both
     being zero.
     The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     of (A,B) satisfies
                      A * v(j) = lambda(j) * B * v(j) .
     The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     of (A,B) satisfies
                      u(j)**H * A  = lambda(j) * u(j)**H * B.
     where u(j)**H is the conjugate-transpose of u(j).
    \endverbatim

 * @param[in] BALANC
          BALANC is CHARACTER*1 \n
          Specifies the balance option to be performed. \n
          = 'N':  do not diagonally scale or permute; \n
          = 'P':  permute only; \n
          = 'S':  scale only; \n
          = 'B':  both permute and scale. \n
          Computed reciprocal condition numbers will be for the
          matrices after permuting and/or balancing. Permuting does
          not change condition numbers (in exact arithmetic), but
          balancing does. \n
 * @param[in] JOBVL
          JOBVL is CHARACTER*1 \n
          = 'N':  do not compute the left generalized eigenvectors; \n
          = 'V':  compute the left generalized eigenvectors. \n
 * @param[in] JOBVR
          JOBVR is CHARACTER*1 \n
          = 'N':  do not compute the right generalized eigenvectors; \n
          = 'V':  compute the right generalized eigenvectors. \n
 * @param[in] SENSE
          SENSE is CHARACTER*1 \n
          Determines which reciprocal condition numbers are computed. \n
          = 'N': none are computed; \n
          = 'E': computed for eigenvalues only; \n
          = 'V': computed for eigenvectors only; \n
          = 'B': computed for eigenvalues and eigenvectors. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, VL, and VR.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the matrix A in the pair (A,B). \n
          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
          or both, then A contains the first part of the real Schur
          form of the "balanced" versions of the input A and B. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the matrix B in the pair (A,B). \n
          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
          or both, then B contains the second part of the real Schur
          form of the "balanced" versions of the input A and B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B.  LDB >= fla_max(1,N). \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
          the j-th eigenvalue is real; if positive, then the j-th and
          (j+1)-st eigenvalues are a complex conjugate pair, with
          ALPHAI(j+1) negative. \n
 \n
          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
          may easily over- or underflow, and BETA(j) may even be zero.
          Thus, the user should avoid naively computing the ratio
          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less
          than and usually comparable with norm(A) in magnitude, and
          BETA always less than and usually comparable with norm(B). \n
 * @param[out] VL
          VL is REAL array, dimension (LDVL,N) \n
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order as
          their eigenvalues. If the j-th eigenvalue is real, then
          u(j) = VL(:,j), the j-th column of VL. If the j-th and
          (j+1)-th eigenvalues form a complex conjugate pair, then
          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
          Each eigenvector will be scaled so the largest component have
          abs(real part) + abs(imag. part) = 1. \n
          Not referenced if JOBVL = 'N'. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the matrix VL. LDVL >= 1, and
          if JOBVL = 'V', LDVL >= N. \n
 * @param[out] VR
          VR is REAL array, dimension (LDVR,N) \n
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order as
          their eigenvalues. If the j-th eigenvalue is real, then
          v(j) = VR(:,j), the j-th column of VR. If the j-th and
          (j+1)-th eigenvalues form a complex conjugate pair, then
          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
          Each eigenvector will be scaled so the largest component have
          abs(real part) + abs(imag. part) = 1. \n
          Not referenced if JOBVR = 'N'. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the matrix VR. LDVR >= 1, and
          if JOBVR = 'V', LDVR >= N. \n
 * @param[out] ILO
          ILO is INTEGER \n
 * @param[out] IHI
          IHI is INTEGER \n
          ILO and IHI are integer values such that on exit
          A(i,j) = 0 and B(i,j) = 0 if i > j and
          j = 1,...,ILO-1 or i = IHI+1,...,N.
          If BALANC = 'N' or 'S', ILO = 1 and IHI = N. \n
 * @param[out] LSCALE
          LSCALE is REAL array, dimension (N) \n
          Details of the permutations and scaling factors applied
          to the left side of A and B.  If PL(j) is the index of the
          row interchanged with row j, and DL(j) is the scaling
          factor applied to row j, then \n
            LSCALE(j) = PL(j)  for j = 1,...,ILO-1 \n
                      = DL(j)  for j = ILO,...,IHI \n
                      = PL(j)  for j = IHI+1,...,N. \n
          The order in which the interchanges are made is N to IHI+1,
          then 1 to ILO-1. \n
 * @param[out] RSCALE
          RSCALE is REAL array, dimension (N) \n
          Details of the permutations and scaling factors applied
          to the right side of A and B.  If PR(j) is the index of the
          column interchanged with column j, and DR(j) is the scaling
          factor applied to column j, then \n
            RSCALE(j) = PR(j)  for j = 1,...,ILO-1 \n
                      = DR(j)  for j = ILO,...,IHI \n
                      = PR(j)  for j = IHI+1,...,N \n
          The order in which the interchanges are made is N to IHI+1,
          then 1 to ILO-1. \n
 * @param[out] ABNRM
          ABNRM is REAL \n
          The one-norm of the balanced matrix A. \n
 * @param[out] BBNRM
          BBNRM is REAL \n
          The one-norm of the balanced matrix B. \n
 * @param[out] RCONDE
          RCONDE is REAL array, dimension (N) \n
          If SENSE = 'E' or 'B', the reciprocal condition numbers of
          the eigenvalues, stored in consecutive elements of the array.
          For a complex conjugate pair of eigenvalues two consecutive
          elements of RCONDE are set to the same value. Thus RCONDE(j),
          RCONDV(j), and the j-th columns of VL and VR all correspond
          to the j-th eigenpair. \n
          If SENSE = 'N' or 'V', RCONDE is not referenced. \n
 * @param[out] RCONDV
          RCONDV is REAL array, dimension (N) \n
          If SENSE = 'V' or 'B', the estimated reciprocal condition
          numbers of the eigenvectors, stored in consecutive elements
          of the array. For a complex eigenvector two consecutive
          elements of RCONDV are set to the same value. If the
          eigenvalues cannot be reordered to compute RCONDV(j),
          RCONDV(j) is set to 0; this can only occur when the true
          value would be very small anyway. \n
          If SENSE = 'N' or 'E', RCONDV is not referenced. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,2*N). \n
          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V',
          LWORK >= fla_max(1,6*N). \n
          If SENSE = 'E', LWORK >= fla_max(1,10*N). \n
          If SENSE = 'V' or 'B', LWORK >= 2*N*N+8*N+16. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N+6) \n
          If SENSE = 'E', IWORK is not referenced. \n
 * @param[out]	BWORK
          BWORK is LOGICAL array, dimension (N) \n
          If SENSE = 'N', BWORK is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1,...,N: \n
                The QZ iteration failed.  No eigenvectors have been
                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
                should be correct for j=INFO+1,...,N. \n
          > N:  =N+1: other than QZ iteration failed in SHGEQZ. \n
                =N+2: error return from STGEVC. \n

 *  * */
    template <typename T>
    void ggevx(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, T *a, integer *lda,
               T *b, integer *ldb, T *alphar, T *alphai, T *beta, T *vl, integer *ldvl, T *vr,
               integer *ldvr, integer *ilo, integer *ihi, T *lscale, T *rscale, T *abnrm, T *bbnrm,
               T *rconde, T *rcondv, T *work, integer *lwork, integer *iwork, logical *bwork,
               integer *info)
    {
        ggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr,
              ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork,
              bwork, info);
    }
    template <typename T, typename Ta>
    void ggevx(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, T *a, integer *lda,
               T *b, integer *ldb, T *alpha, T *beta, T *vl, integer *ldvl, T *vr, integer *ldvr,
               integer *ilo, integer *ihi, Ta *lscale, Ta *rscale, Ta *abnrm, Ta *bbnrm, Ta *rconde,
               Ta *rcondv, T *work, integer *lwork, Ta *rwork, integer *iwork, logical *bwork,
               integer *info)
    {
        ggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo,
              ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork,
              info);
    }
    /** @}*/ // end of ggevx

    /** @defgroup gges3 gges3
     * @ingroup NYS
     * @{
     */
    /*! @brief GGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
 vectors for GE matrices (blocked algorithm)
 * @details
 * \b Purpose:
    \verbatim
     GGES3 computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     the generalized eigenvalues, the generalized real Schur form (S,T),
     optionally, the left and/or right matrices of Schur vectors (VSL and
     VSR). This gives the generalized Schur factorization

              (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T)

     Optionally, it also orders the eigenvalues so that a selected cluster
     of eigenvalues appears in the leading diagonal blocks of the upper
     quasi-triangular matrix S and the upper triangular matrix T.The
     leading columns of VSL and VSR then form an orthonormal basis for the
     corresponding left and right eigenspaces (deflating subspaces).

     (If only the generalized eigenvalues are needed, use the driver
     SGGEV instead, which is faster.)

     A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     usually represented as the pair (alpha,beta), as there is a
     reasonable interpretation for beta=0 or both being zero.

     A pair of matrices (S,T) is in generalized real Schur form if T is
     upper triangular with non-negative diagonal and S is block upper
     triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     to real generalized eigenvalues, while 2-by-2 blocks of S will be
     "standardized" by making the corresponding elements of T have the
     form:
             [  a  0  ]
             [  0  b  ]

     and the pair of corresponding 2-by-2 blocks in S and T will have a
     complex conjugate pair of generalized eigenvalues.
    \endverbatim

 * @param[in] JOBVSL
          JOBVSL is CHARACTER*1 \n
          = 'N':  do not compute the left Schur vectors; \n
          = 'V':  compute the left Schur vectors. \n
 * @param[in] JOBVSR
          JOBVSR is CHARACTER*1 \n
          = 'N':  do not compute the right Schur vectors; \n
          = 'V':  compute the right Schur vectors. \n
 * @param[in] SORT
          SORT is CHARACTER*1 \n
          Specifies whether or not to order the eigenvalues on the
          diagonal of the generalized Schur form. \n
          = 'N':  Eigenvalues are not ordered; \n
          = 'S':  Eigenvalues are ordered (see SELCTG); \n
 * @param[in] SELCTG
          SELCTG is a LOGICAL FUNCTION of three REAL arguments \n
          SELCTG must be declared EXTERNAL in the calling subroutine.
          If SORT = 'N', SELCTG is not referenced.
          If SORT = 'S', SELCTG is used to select eigenvalues to sort
          to the top left of the Schur form. \n
          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
          one of a complex conjugate pair of eigenvalues is selected,
          then both complex eigenvalues are selected. \n
 \n
          Note that in the ill-conditioned case, a selected complex
          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),
          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2
          in this case. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, VSL, and VSR.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the first of the pair of matrices.
          On exit, A has been overwritten by its generalized Schur
          form S. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the second of the pair of matrices.
          On exit, B has been overwritten by its generalized Schur
          form T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B.  LDB >= fla_max(1,N). \n
 * @param[out] SDIM
          SDIM is INTEGER \n
          If SORT = 'N', SDIM = 0. \n
          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
          for which SELCTG is true.  (Complex conjugate pairs for which
          SELCTG is true for either eigenvalue count as 2.) \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
          and  BETA(j),j=1,...,N are the diagonals of the complex Schur
          form (S,T) that would result if the 2-by-2 diagonal blocks of
          the real Schur form of (A,B) were further reduced to
          triangular form using 2-by-2 complex unitary transformations.
          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
          positive, then the j-th and (j+1)-st eigenvalues are a
          complex conjugate pair, with ALPHAI(j+1) negative. \n
 \n
          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
          may easily over- or underflow, and BETA(j) may even be zero.
          Thus, the user should avoid naively computing the ratio.
          However, ALPHAR and ALPHAI will be always less than and
          usually comparable with norm(A) in magnitude, and BETA always
          less than and usually comparable with norm(B). \n
 * @param[out] VSL
          VSL is REAL array, dimension (LDVSL,N) \n
          If JOBVSL = 'V', VSL will contain the left Schur vectors.
          Not referenced if JOBVSL = 'N'. \n
 * @param[in] LDVSL
          LDVSL is INTEGER \n
          The leading dimension of the matrix VSL. LDVSL >=1, and
          if JOBVSL = 'V', LDVSL >= N. \n
 * @param[out] VSR
          VSR is REAL array, dimension (LDVSR,N) \n
          If JOBVSR = 'V', VSR will contain the right Schur vectors.
          Not referenced if JOBVSR = 'N'. \n
 * @param[in] LDVSR
          LDVSR is INTEGER \n
          The leading dimension of the matrix VSR. LDVSR >= 1, and
          if JOBVSR = 'V', LDVSR >= N. \n
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
 * @param[out]	BWORK
          BWORK is LOGICAL array, dimension (N) \n
          Not referenced if SORT = 'N'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1,...,N: \n
                The QZ iteration failed.  (A,B) are not in Schur
                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
                be correct for j=INFO+1,...,N. \n
          > N:  =N+1: other than QZ iteration failed in SHGEQZ. \n
                =N+2: after reordering, roundoff changed values of
                      some complex eigenvalues so that leading
                      eigenvalues in the Generalized Schur form no
                      longer satisfy SELCTG=.TRUE.  This could also
                      be caused due to scaling. \n
                =N+3: reordering failed in STGSEN. \n

 *  * */
    template <typename T>
    void gges3(char *jobvsl, char *jobvsr, char *sort, void *selctg, integer *n, T *a, integer *lda,
               T *b, integer *ldb, integer *sdim, T *alphar, T *alphai, T *beta, T *vsl,
               integer *ldvsl, T *vsr, integer *ldvsr, T *work, integer *lwork, logical *bwork,
               integer *info)
    {
        gges3(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl,
              ldvsl, vsr, ldvsr, work, lwork, bwork, info);
    }
    template <typename T, typename Ta>
    void gges3(char *jobvsl, char *jobvsr, char *sort, void *selctg, integer *n, T *a, integer *lda,
               T *b, integer *ldb, integer *sdim, T *alpha, T *beta, T *vsl, integer *ldvsl, T *vsr,
               integer *ldvsr, T *work, integer *lwork, Ta *rwork, logical *bwork, integer *info)
    {
        gges3(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr,
              ldvsr, work, lwork, rwork, bwork, info);
    }
    /** @}*/ // end of gges3

    /** @defgroup gges gges
     * @ingroup NYS
     * @{
     */
    /*! @brief GGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
 vectors for GE matrices
 * @details
 * \b Purpose:
    \verbatim
     GGES computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     the generalized eigenvalues, the generalized real Schur form (S,T),
     optionally, the left and/or right matrices of Schur vectors (VSL and
     VSR). This gives the generalized Schur factorization

              (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T)

     Optionally, it also orders the eigenvalues so that a selected cluster
     of eigenvalues appears in the leading diagonal blocks of the upper
     quasi-triangular matrix S and the upper triangular matrix T.The
     leading columns of VSL and VSR then form an orthonormal basis for the
     corresponding left and right eigenspaces (deflating subspaces).

     (If only the generalized eigenvalues are needed, use the driver
     SGGEV instead, which is faster.)

     A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     usually represented as the pair (alpha,beta), as there is a
     reasonable interpretation for beta=0 or both being zero.

     A pair of matrices (S,T) is in generalized real Schur form if T is
     upper triangular with non-negative diagonal and S is block upper
     triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     to real generalized eigenvalues, while 2-by-2 blocks of S will be
     "standardized" by making the corresponding elements of T have the
     form:
             [  a  0  ]
             [  0  b  ]

     and the pair of corresponding 2-by-2 blocks in S and T will have a
     complex conjugate pair of generalized eigenvalues.
    \endverbatim

 * @param[in] JOBVSL
          JOBVSL is CHARACTER*1 \n
          = 'N':  do not compute the left Schur vectors; \n
          = 'V':  compute the left Schur vectors. \n
 * @param[in] JOBVSR
          JOBVSR is CHARACTER*1 \n
          = 'N':  do not compute the right Schur vectors; \n
          = 'V':  compute the right Schur vectors. \n
 * @param[in] SORT
          SORT is CHARACTER*1 \n
          Specifies whether or not to order the eigenvalues on the
          diagonal of the generalized Schur form. \n
          = 'N':  Eigenvalues are not ordered; \n
          = 'S':  Eigenvalues are ordered (see SELCTG); \n
 * @param[in] SELCTG
          SELCTG is a LOGICAL FUNCTION of three REAL arguments
          SELCTG must be declared EXTERNAL in the calling subroutine. \n
          If SORT = 'N', SELCTG is not referenced. \n
          If SORT = 'S', SELCTG is used to select eigenvalues to sort
          to the top left of the Schur form. \n
          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
          one of a complex conjugate pair of eigenvalues is selected,
          then both complex eigenvalues are selected. \n
 \n
          Note that in the ill-conditioned case, a selected complex
          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),
          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2
          in this case. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, VSL, and VSR.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the first of the pair of matrices. \n
          On exit, A has been overwritten by its generalized Schur
          form S. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the second of the pair of matrices. \n
          On exit, B has been overwritten by its generalized Schur
          form T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B.  LDB >= fla_max(1,N). \n
 * @param[out] SDIM
          SDIM is INTEGER \n
          If SORT = 'N', SDIM = 0. \n
          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
          for which SELCTG is true.  (Complex conjugate pairs for which
          SELCTG is true for either eigenvalue count as 2.) \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
          and  BETA(j),j=1,...,N are the diagonals of the complex Schur
          form (S,T) that would result if the 2-by-2 diagonal blocks of
          the real Schur form of (A,B) were further reduced to
          triangular form using 2-by-2 complex unitary transformations.
          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
          positive, then the j-th and (j+1)-st eigenvalues are a
          complex conjugate pair, with ALPHAI(j+1) negative. \n
 \n
          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
          may easily over- or underflow, and BETA(j) may even be zero.
          Thus, the user should avoid naively computing the ratio.
          However, ALPHAR and ALPHAI will be always less than and
          usually comparable with norm(A) in magnitude, and BETA always
          less than and usually comparable with norm(B). \n
 * @param[out] VSL
          VSL is REAL array, dimension (LDVSL,N) \n
          If JOBVSL = 'V', VSL will contain the left Schur vectors.
          Not referenced if JOBVSL = 'N'. \n
 * @param[in] LDVSL
          LDVSL is INTEGER \n
          The leading dimension of the matrix VSL. LDVSL >=1, and
          if JOBVSL = 'V', LDVSL >= N. \n
 * @param[out] VSR
          VSR is REAL array, dimension (LDVSR,N) \n
          If JOBVSR = 'V', VSR will contain the right Schur vectors.
          Not referenced if JOBVSR = 'N'. \n
 * @param[in] LDVSR
          LDVSR is INTEGER \n
          The leading dimension of the matrix VSR. LDVSR >= 1, and
          if JOBVSR = 'V', LDVSR >= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If N = 0, LWORK >= 1, else LWORK >= fla_max(8*N,6*N+16).
          For good performance , LWORK must generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	BWORK
          BWORK is LOGICAL array, dimension (N) \n
          Not referenced if SORT = 'N'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1,...,N:
                The QZ iteration failed.  (A,B) are not in Schur
                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
                be correct for j=INFO+1,...,N. \n
          > N:  =N+1: other than QZ iteration failed in SHGEQZ. \n
                =N+2: after reordering, roundoff changed values of
                      some complex eigenvalues so that leading
                      eigenvalues in the Generalized Schur form no
                      longer satisfy SELCTG=.TRUE.  This could also
                      be caused due to scaling. \n
                =N+3: reordering failed in STGSEN. \n

 *  * */
    template <typename T>
    void gges(char *jobvsl, char *jobvsr, char *sort, void *selctg, integer *n, T *a, integer *lda,
              T *b, integer *ldb, integer *sdim, T *alphar, T *alphai, T *beta, T *vsl,
              integer *ldvsl, T *vsr, integer *ldvsr, T *work, integer *lwork, logical *bwork,
              integer *info)
    {
        gges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl,
             ldvsl, vsr, ldvsr, work, lwork, bwork, info);
    }
    template <typename T, typename Ta>
    void gges(char *jobvsl, char *jobvsr, char *sort, void *selctg, integer *n, T *a, integer *lda,
              T *b, integer *ldb, integer *sdim, T *alpha, T *beta, T *vsl, integer *ldvsl, T *vsr,
              integer *ldvsr, T *work, integer *lwork, Ta *rwork, logical *bwork, integer *info)
    {
        gges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr,
             ldvsr, work, lwork, rwork, bwork, info);
    }
    /** @}*/ // end of gges

    /** @defgroup ggesx ggesx
     * @ingroup NYS
     * @{
     */
    /*! @brief GGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
 vectors for GE matrices

 * @details
 * \b Purpose:
    \verbatim
      GGESX computes for a pair of N-by-N real nonsymmetric matrices
     (A,B), the generalized eigenvalues, the real Schur form (S,T), and,
     optionally, the left and/or right matrices of Schur vectors (VSL and
     VSR).  This gives the generalized Schur factorization

          (A,B) = ( (VSL) S (VSR)**T, (VSL) T (VSR)**T)

     Optionally, it also orders the eigenvalues so that a selected cluster
     of eigenvalues appears in the leading diagonal blocks of the upper
     quasi-triangular matrix S and the upper triangular matrix T; computes
     a reciprocal condition number for the average of the selected
     eigenvalues (RCONDE); and computes a reciprocal condition number for
     the right and left deflating subspaces corresponding to the selected
     eigenvalues (RCONDV). The leading columns of VSL and VSR then form
     an orthonormal basis for the corresponding left and right eigenspaces
     (deflating subspaces).

     A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     usually represented as the pair (alpha,beta), as there is a
     reasonable interpretation for beta=0 or for both being zero.

     A pair of matrices (S,T) is in generalized real Schur form if T is
     upper triangular with non-negative diagonal and S is block upper
     triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     to real generalized eigenvalues, while 2-by-2 blocks of S will be
     "standardized" by making the corresponding elements of T have the
     form:
             [  a  0  ]
             [  0  b  ]

     and the pair of corresponding 2-by-2 blocks in S and T will have a
     complex conjugate pair of generalized eigenvalues.
    \endverbatim

 * @param[in] JOBVSL
          JOBVSL is CHARACTER*1 \n
          = 'N':  do not compute the left Schur vectors; \n
          = 'V':  compute the left Schur vectors. \n
 * @param[in] JOBVSR
          JOBVSR is CHARACTER*1 \n
          = 'N':  do not compute the right Schur vectors; \n
          = 'V':  compute the right Schur vectors. \n
 * @param[in] SORT
          SORT is CHARACTER*1 \n
          Specifies whether or not to order the eigenvalues on the
          diagonal of the generalized Schur form. \n
          = 'N':  Eigenvalues are not ordered; \n
          = 'S':  Eigenvalues are ordered (see SELCTG). \n
 * @param[in] SELCTG
          SELCTG is a LOGICAL FUNCTION of three REAL arguments \n
          SELCTG must be declared EXTERNAL in the calling subroutine. \n
          If SORT = 'N', SELCTG is not referenced. \n
          If SORT = 'S', SELCTG is used to select eigenvalues to sort
          to the top left of the Schur form. \n
          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
          one of a complex conjugate pair of eigenvalues is selected,
          then both complex eigenvalues are selected. \n
          Note that a selected complex eigenvalue may no longer satisfy
          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) = .TRUE. after ordering,
          since ordering may change the value of complex eigenvalues
          (especially if the eigenvalue is ill-conditioned), in this
          case INFO is set to N+3. \n
 * @param[in] SENSE
          SENSE is CHARACTER*1 \n
          Determines which reciprocal condition numbers are computed. \n
          = 'N':  None are computed; \n
          = 'E':  Computed for average of selected eigenvalues only; \n
          = 'V':  Computed for selected deflating subspaces only; \n
          = 'B':  Computed for both. \n
          If SENSE = 'E', 'V', or 'B', SORT must equal 'S'. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, VSL, and VSR.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the first of the pair of matrices. \n
          On exit, A has been overwritten by its generalized Schur
          form S. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the second of the pair of matrices. \n
          On exit, B has been overwritten by its generalized Schur
          form T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B.  LDB >= fla_max(1,N). \n
 * @param[out] SDIM
          SDIM is INTEGER \n
          If SORT = 'N', SDIM = 0. \n
          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
          for which SELCTG is true.  (Complex conjugate pairs for which
          SELCTG is true for either eigenvalue count as 2.) \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i
          and BETA(j),j=1,...,N  are the diagonals of the complex Schur
          form (S,T) that would result if the 2-by-2 diagonal blocks of
          the real Schur form of (A,B) were further reduced to
          triangular form using 2-by-2 complex unitary transformations.
          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
          positive, then the j-th and (j+1)-st eigenvalues are a
          complex conjugate pair, with ALPHAI(j+1) negative. \n
 \n
          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
          may easily over- or underflow, and BETA(j) may even be zero.
          Thus, the user should avoid naively computing the ratio.
          However, ALPHAR and ALPHAI will be always less than and
          usually comparable with norm(A) in magnitude, and BETA always
          less than and usually comparable with norm(B). \n
 * @param[out] VSL
          VSL is REAL array, dimension (LDVSL,N) \n
          If JOBVSL = 'V', VSL will contain the left Schur vectors.
          Not referenced if JOBVSL = 'N'. \n
 * @param[in] LDVSL
          LDVSL is INTEGER \n
          The leading dimension of the matrix VSL. LDVSL >=1, and
          if JOBVSL = 'V', LDVSL >= N. \n
 * @param[out] VSR
          VSR is REAL array, dimension (LDVSR,N) \n
          If JOBVSR = 'V', VSR will contain the right Schur vectors.
          Not referenced if JOBVSR = 'N'. \n
 * @param[in] LDVSR
          LDVSR is INTEGER \n
          The leading dimension of the matrix VSR. LDVSR >= 1, and
          if JOBVSR = 'V', LDVSR >= N. \n
 * @param[out] RCONDE
          RCONDE is REAL array, dimension ( 2) \n
          If SENSE = 'E' or 'B', RCONDE(1) and RCONDE(2) contain the
          reciprocal condition numbers for the average of the selected
          eigenvalues.
          Not referenced if SENSE = 'N' or 'V'. \n
 * @param[out] RCONDV
          RCONDV is REAL array, dimension ( 2) \n
          If SENSE = 'V' or 'B', RCONDV(1) and RCONDV(2) contain the
          reciprocal condition numbers for the selected deflating
          subspaces.
          Not referenced if SENSE = 'N' or 'E'. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If N = 0, LWORK >= 1, else if SENSE = 'E', 'V', or 'B',
          LWORK >= fla_max( 8*N, 6*N+16, 2*SDIM*(N-SDIM) ), else
          LWORK >= fla_max( 8*N, 6*N+16 ). \n
          Note that 2*SDIM*(N-SDIM) <= N*N/2.
          Note also that an error is only returned if
          LWORK < fla_max( 8*N, 6*N+16), but if SENSE = 'E' or 'V' or 'B'
          this may not be large enough. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the bound on the optimal size of the WORK
          array and the minimum size of the IWORK array, returns these
          values as the first entries of the WORK and IWORK arrays, and
          no error message related to LWORK or LIWORK is issued by
          XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          If SENSE = 'N' or N = 0, LIWORK >= 1, otherwise
          LIWORK >= N+6. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the bound on the optimal size of the
          WORK array and the minimum size of the IWORK array, returns
          these values as the first entries of the WORK and IWORK
          arrays, and no error message related to LWORK or LIWORK is
          issued by XERBLA. \n
 * @param[out]	BWORK
          BWORK is LOGICAL array, dimension (N) \n
          Not referenced if SORT = 'N'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1,...,N: \n
                The QZ iteration failed.  (A,B) are not in Schur
                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
                be correct for j=INFO+1,...,N. \n
          > N:  =N+1: other than QZ iteration failed in SHGEQZ \n
                =N+2: after reordering, roundoff changed values of
                      some complex eigenvalues so that leading
                      eigenvalues in the Generalized Schur form no
                      longer satisfy SELCTG=.TRUE.  This could also
                      be caused due to scaling. \n
                =N+3: reordering failed in STGSEN. \n

 *  * */
    template <typename T>
    void ggesx(char *jobvsl, char *jobvsr, char *sort, void *selctg, char *sense, integer *n, T *a,
               integer *lda, T *b, integer *ldb, integer *sdim, T *alphar, T *alphai, T *beta,
               T *vsl, integer *ldvsl, T *vsr, integer *ldvsr, T *rconde, T *rcondv, T *work,
               integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info)
    {
        ggesx(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta,
              vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
    }
    template <typename T, typename Ta>
    void ggesx(char *jobvsl, char *jobvsr, char *sort, void *selctg, char *sense, integer *n, T *a,
               integer *lda, T *b, integer *ldb, integer *sdim, T *alpha, T *beta, T *vsl,
               integer *ldvsl, T *vsr, integer *ldvsr, Ta *rconde, Ta *rcondv, T *work,
               integer *lwork, Ta *rwork, integer *iwork, integer *liwork, logical *bwork,
               integer *info)
    {
        ggesx(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl,
              vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, info);
    }
    /** @}*/ // end of ggesx

    /** @defgroup gebal gebal
     * @ingroup NYS
     * @{
     */
    /*! @brief \b GEBAL It balances a real matrix
 *
 *  @details
 *  \b Purpose:
    \verbatim
    GEBAL balances a general real matrix A.  This involves, first,
    permuting A by a similarity transformation to isolate eigenvalues
    in the first 1 to ILO-1 and last IHI+1 to N elements on the
    diagonal; and second, applying a diagonal similarity transformation
    to rows and columns ILO to IHI to make the rows and columns as
    close in norm as possible.  Both steps are optional.

    Balancing may reduce the 1-norm of the matrix, and improve the
    accuracy of the computed eigenvalues and/or eigenvectors.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          Specifies the operations to be performed on A: \n
          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
                  for i = 1,...,N; \n
          = 'P':  permute only; \n
          = 'S':  scale only; \n
          = 'B':  both permute and scale. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the input matrix A. \n
          On exit,  A is overwritten by the balanced matrix.
          If JOB = 'N', A is not referenced.
          See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] ILO
          ILO is INTEGER \n
 * @param[out] IHI
          IHI is INTEGER \n
          ILO and IHI are set to integers such that on exit
          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
          If JOB = 'N' or 'S', ILO = 1 and IHI = N. \n
 * @param[out] SCALE
          SCALE is REAL array, dimension (N) \n
          Details of the permutations and scaling factors applied to
          A.  If P(j) is the index of the row and column interchanged
          with row and column j and D(j) is the scaling factor
          applied to row and column j, then \n
          SCALE(j) = P(j)    for j = 1,...,ILO-1 \n
                   = D(j)    for j = ILO,...,IHI \n
                   = P(j)    for j = IHI+1,...,N. \n
          The order in which the interchanges are made is N to IHI+1, \n
          then 1 to ILO-1. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void gebal(char *job, integer *n, T *a, integer *lda, integer *ilo, integer *ihi, T *scale,
               integer *info)
    {
        gebal(job, n, a, lda, ilo, ihi, scale, info);
    }
    template <typename T, typename Ta>
    void gebal(char *job, integer *n, T *a, integer *lda, integer *ilo, integer *ihi, Ta *scale,
               integer *info)
    {
        gebal(job, n, a, lda, ilo, ihi, scale, info);
    }
    /** @}*/ // end of gebal

    /** @defgroup gehrd gehrd
     * @ingroup NYS
     * @{
     */
    /*! @brief Reduction to upper Hessenberg form
    *
    * @details
    * \b Purpose:
    * \verbatim
    Reduction of a real general matrix a to upper Hessenberg form H by an orthogonal
    similarity transformation:
    Q**T * A * Q = H .
    \endverbatim

    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in] ilo
              ilo is integer* \n
    * @param[in] ihi
              ihi is integer* \n
              It is assumed that A is already upper triangular in rows
              and columns 1:ilo-1 and ihi+1:n. \n
              ilo and ihi are normally set by a previous call to SGEBAL;
              otherwise they should be set to 1 and N respectively. See Further Details. \n
              1 <= ilo <= ihi <= N, if N > 0; ilo=1 and ihi=0, if N=0. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the n-by-n general matrix to be reduced. \n
              On exit, the upper triangle and the first subdiagonal of A
              are overwritten with the upper Hessenberg matrix H, and the
              elements below the first subdiagonal, with the array tau,
              represent the orthogonal matrix Q as a product of elementary
              reflectors. See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] tau
              tau is float/double/COMPLEX/COMPLEX*16 array, dimension (n-1) \n
              The scalar factors of the elementary reflectors (see Further
              Details). Elements 1:ilo-1 and ihi:N-1 of tau are set to
              zero.
              *
              * \n
              * **Further Details**
              * \verbatim
                      The matrix Q is represented as a product of (ihi-ilo) elementary reflectors

                      Q = H(ilo) H(ilo+1) . . . H(ihi-1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**T

                      where tau is a real scalar, and V is a real vector with V(1:i) = 0, V(i+1) = 1
 and V(ihi+1:n) = 0; V(i+2:ihi) is stored on exit in A(i+2:ihi,i), and tau in tau(i).

                      The contents of A are illustrated by the following example, with n = 7, ilo =
 2 and ihi = 6: on entry,   on exit,

                      ( a   a   a   a   a   a   a) (  a   a   h   h   h   h   a)
                      (  a   a   a   a   a   a)    (   a   h   h   h   h   a)
                      (  a   a   a   a   a   a)    (   h   h   h   h   h   h)
                      (  a   a   a   a   a   a)    (   v2  h   h   h   h   h)
                      (  a   a   a   a   a   a)    (   v2  v3  h   h   h   h)
                      (  a   a   a   a   a   a)    (   v2  v3  v4  h   h   h)
                      ( a) (  a)
                      where,
                      a denotes an element of the original matrix a,
                      h denotes a modified element of the upper Hessenberg matrix H,
                      vi denotes an element of the vector defining H(i).
              \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (LWORK) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The length of the array WORK.  LWORK >= fla_max(1,N). \n
              For good performance, LWORK should generally be larger. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void gehrd(integer *n, integer *ilo, integer *ihi, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        gehrd(n, ilo, ihi, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of gehrd

    /** @defgroup gehd2 gehd2
     * @ingroup NYS
     * @{
     */
    /*! @brief Reduction to upper Hessenberg form using an unblocked algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real general matrix a to upper Hessenberg form H by an orthogonal
        similarity transformation:
        Q**T * A * Q = H .
    \endverbatim

    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in] ilo
              ilo is integer* \n
    * @param[in] ihi
              ihi is integer* \n
              It is assumed that A is already upper triangular in rows
              and columns 1:ilo-1 and ihi+1:n. \n
              ilo and ihi are normally set by a previous call to SGEBAL;
              otherwise they should be set to 1 and N respectively. See Further Details. \n
              1 <= ilo <= ihi <= fla_max(1,n). \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the n-by-n general matrix to be reduced. \n
              On exit, the upper triangle and the first subdiagonal of A
              are overwritten with the upper Hessenberg matrix H, and the
              elements below the first subdiagonal, with the array tau,
              represent the orthogonal matrix Q as a product of elementary
              reflectors. See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] tau
              tau is float/double/COMPLEX/COMPLEX*16 array, dimension (n-1) \n
              The scalar factors of the elementary reflectors (see Further
              Details).
              *
              * \n
              * **Further Details**
              * \verbatim
                      The matrix Q is represented as a product of (ihi-ilo) elementary reflectors

                      Q = H(ilo) H(ilo+1) . . . H(ihi-1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**T

                      where tau is a real scalar, and V is a real vector with V(1:i) = 0, V(i+1) = 1
    and V(ihi+1:n) = 0; V(i+2:ihi) is stored on exit in A(i+2:ihi,i), and tau in tau(i).

                      The contents of A are illustrated by the following example, with n = 7, ilo =
    2 and ihi = 6: on entry,   on exit,

                      ( a   a   a   a   a   a   a) (  a   a   h   h   h   h   a)
                      (  a   a   a   a   a   a) (   a   h   h   h   h   a)
                      (  a   a   a   a   a   a) (   h   h   h   h   h   h)
                      (  a   a   a   a   a   a) (   v2  h   h   h   h   h)
                      (  a   a   a   a   a   a) (   v2  v3  h   h   h   h)
                      (  a   a   a   a   a   a) (   v2  v3  v4  h   h   h)
                      ( a) (  a)
                      where,
                      a denotes an element of the original matrix a,
                      h denotes a modified element of the upper Hessenberg matrix H,
                      vi denotes an element of the vector defining H(i).
              \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (N) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void gehd2(integer *n, integer *ilo, integer *ihi, T *a, integer *lda, T *tau, T *work,
               integer *info)
    {
        gehd2(n, ilo, ihi, a, lda, tau, work, info);
    }
    /** @}*/ // end of gehd2

    /** @defgroup lahr2 lahr2
     * @ingroup NYS
     * @{
     */
    /*! @brief LAHR2 reduces the specified number of first columns of a general \n
     rectangular matrix A so that elements below the specified subdiagonal  \n
     are zero, and   returns auxiliary matrices which are needed to apply the \n
     transformation to the unreduced part of A
 * @details
 * \b Purpose:
    \verbatim
    LAHR2 reduces the first NB columns of A real general n-BY-(n-k+1)
    matrix A so that elements below the k-th subdiagonal are zero. The
    reduction is performed by an orthogonal similarity transformation
    Q**T * A * Q. The routine   returns the matrices V and T which determine
    Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.

    This is an auxiliary routine called by GEHRD.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. \n
 * @param[in] K
          K is INTEGER \n
          The offset for the reduction. Elements below the k-th
          subdiagonal in the first NB columns are reduced to zero.
          K < N. \n
 * @param[in] NB
          NB is INTEGER \n
          The number of columns to be reduced. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N-K+1) \n
          On entry, the n-by-(n-k+1) general matrix A. \n
          On exit, the elements on and above the k-th subdiagonal in
          the first NB columns are overwritten with the corresponding
          elements of the reduced matrix; the elements below the k-th
          subdiagonal, with the array TAU, represent the matrix Q as a
          product of elementary reflectors. The other columns of A are
          unchanged. See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] TAU
          TAU is REAL array, dimension (NB) \n
          The scalar factors of the elementary reflectors. See Further
          Details. \n
 * @param[out] T
          T is REAL array, dimension (LDT,NB) \n
          The upper triangular matrix T. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= NB. \n
 * @param[out] Y
          Y is REAL array, dimension (LDY,NB) \n
          The n-by-nb matrix Y. \n
 * @param[in] LDY
          LDY is INTEGER \n
          The leading dimension of the array Y. LDY >= N.  \n

 * */
    template <typename T>
    void lahr2(integer *n, integer *k, integer *nb, T *a, integer *lda, T *tau, T *t, integer *ldt,
               T *y, integer *ldy)
    {
        lahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy);
    }
    /** @}*/ // end of lahr2

    /** @defgroup ghr {un,or}ghr
     * @ingroup NYS
     * @{
     */
    /*! @brief ORGHR generates a real orthogonal matrix Q

 * @details
 * \b Purpose:
    \verbatim
     ORGHR generates a real orthogonal matrix Q which is defined as the
     product of IHI-ILO elementary reflectors of order N, as returned by
     SGEHRD:

     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix Q. N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          ILO and IHI must have the same values as in the previous call
          of SGEHRD. Q is equal to the unit matrix except in the
          submatrix Q(ilo+1:ihi,ilo+1:ihi). \n
          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the vectors which define the elementary reflectors,
          as returned by SGEHRD. \n
          On exit, the N-by-N orthogonal matrix Q. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[in] TAU
          TAU is REAL array, dimension (N-1) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by SGEHRD. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= IHI-ILO.
          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
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
    template <typename T>
    void orghr(integer *n, integer *ilo, integer *ihi, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        orghr(n, ilo, ihi, a, lda, tau, work, lwork, info);
    }
    template <typename T>
    void unghr(integer *n, integer *ilo, integer *ihi, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        unghr(n, ilo, ihi, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of ghr

    /** @defgroup mhr {un, or}mhr
     * @ingroup NYS
     * @{
     */
    /*! @brief ORMHR overwrites the general real M-by-N matrix C

 * @details
 * \b Purpose:
    \verbatim
    ORMHR overwrites the general real M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'T':      Q**T * C       C * Q**T

    where Q is a real orthogonal matrix of order nq, with nq = m if
    SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
    IHI-ILO elementary reflectors, as returned by SGEHRD:

    Q = H(ilo) H(ilo+1) . . . H(ihi-1).
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          ILO and IHI must have the same values as in the previous call
          of SGEHRD. Q is equal to the unit matrix except in the
          submatrix Q(ilo+1:ihi,ilo+1:ihi). \n
          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and
          ILO = 1 and IHI = 0, if M = 0; \n
          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and
          ILO = 1 and IHI = 0, if N = 0. \n
 * @param[in] A
          A is REAL array, dimension \n
                               (LDA,M) if SIDE = 'L' \n
                               (LDA,N) if SIDE = 'R' \n
          The vectors which define the elementary reflectors, as
          returned by SGEHRD. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. \n
          LDA >= fla_max(1,M) if SIDE = 'L'; LDA >= fla_max(1,N) if SIDE = 'R'. \n
 * @param[in] TAU
          TAU is REAL array, dimension \n
                               (M-1) if SIDE = 'L' \n
                               (N-1) if SIDE = 'R' \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by SGEHRD. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C. \n
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
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

 *  * */
    template <typename T>
    void ormhr(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info);
    }
    template <typename T>
    void unmhr(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @}*/ // end of mhr

    /** @defgroup gebak gebak
     * @ingroup NYS
     * @{
     */
    /*! @brief \b GEBAK Forms right left vectors of real general matrix
 *
 * @details
    \b Purpose:
    *  =============
    \verbatim
    GEBAK forms the right or left eigenvectors of a real general matrix
    by backward transformation on the computed eigenvectors of the
    balanced matrix output by GEBAL.
    \endverbatim

 * @param[in] JOB
        JOB is CHARACTER*1 \n
        Specifies the type of backward transformation required: \n
        = 'N': do nothing, return immediately; \n
        = 'P': do backward transformation for permutation only; \n
        = 'S': do backward transformation for scaling only; \n
        = 'B': do backward transformations for both permutation and
               scaling. \n
        JOB must be the same as the argument JOB supplied to SGEBAL. \n
 * @param[in] SIDE
        SIDE is CHARACTER*1 \n
        = 'R':  V contains right eigenvectors; \n
        = 'L':  V contains left eigenvectors. \n
 * @param[in] N
        N is INTEGER \n
        The number of rows of the matrix V.  N >= 0. \n
 * @param[in] ILO
        ILO is INTEGER \n
 * @param[in] IHI
        IHI is INTEGER \n
        The integers ILO and IHI determined by SGEBAL. \n
        1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. \n
 * @param[in] SCALE
        SCALE is REAL array, dimension (N) \n
        Details of the permutation and scaling factors, as returned
        by SGEBAL. \n
 * @param[in] M
        M is INTEGER \n
        The number of columns of the matrix V.  M >= 0. \n
 * @param[in,out] V
        V is REAL array, dimension (LDV,M) \n
        On entry, the matrix of right or left eigenvectors to be
        transformed, as returned by SHSEIN or STREVC. \n
        On exit, V is overwritten by the transformed eigenvectors. \n
 * @param[in] LDV
        LDV is INTEGER \n
        The leading dimension of the array V. LDV >= fla_max(1,N). \n
 * @param[out]	INFO
        INFO is INTEGER \n
        = 0:  successful exit \n
        < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void gebak(char *job, char *side, integer *n, integer *ilo, integer *ihi, T *scale, integer *m,
               T *v, integer *ldv, integer *info)
    {
        gebak(job, side, n, ilo, ihi, scale, m, v, ldv, info);
    }
    template <typename T, typename Ta>
    void gebak(char *job, char *side, integer *n, integer *ilo, integer *ihi, Ta *scale, integer *m,
               T *v, integer *ldv, integer *info)
    {
        gebak(job, side, n, ilo, ihi, scale, m, v, ldv, info);
    }
    /** @}*/ // end of gebak

    /** @defgroup hseqr hseqr
     * @ingroup NYS
     * @{
     */
    /*! @brief HSEQR computes the eigenvalues of a Hessenberg matrix H

 * @details
 * \b Purpose:
    \verbatim
    HSEQR computes the eigenvalues of a Hessenberg matrix H
    and, optionally, the matrices T and Z from the Schur decomposition
    H = Z T Z**T, where T is an upper quasi-triangular matrix (the
    Schur form), and Z is the orthogonal matrix of Schur vectors.

    Optionally Z may be postmultiplied into an input orthogonal
    matrix Q so that this routine can give the Schur factorization
    of a matrix A which has been reduced to the Hessenberg form H
    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
           = 'E':  compute eigenvalues only; \n
           = 'S':  compute eigenvalues and the Schur form T. \n
 * @param[in] COMPZ
          COMPZ is CHARACTER*1 \n
           = 'N':  no Schur vectors are computed; \n
           = 'I':  Z is initialized to the unit matrix and the matrix Z
                   of Schur vectors of H is returned; \n
           = 'V':  Z must contain an orthogonal matrix Q on entry, and
                   the product Q*Z is returned. \n
 * @param[in] N
          N is INTEGER \n
           The order of the matrix H.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          It is assumed that H is already upper triangular in rows
          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
          set by a previous call to SGEBAL, and then passed to ZGEHRD
          when the matrix output by SGEBAL is reduced to Hessenberg
          form. Otherwise ILO and IHI should be set to 1 and N
          respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
          If N = 0, then ILO = 1 and IHI = 0. \n
 * @param[in,out] H
          H is REAL array, dimension (LDH,N) \n
          On entry, the upper Hessenberg matrix H. \n
          On exit, if INFO = 0 and JOB = 'S', then H contains the
          upper quasi-triangular matrix T from the Schur decomposition
          (the Schur form); 2-by-2 diagonal blocks (corresponding to
          complex conjugate pairs of eigenvalues) are returned in
          standard form, with H(i,i) = H(i+1,i+1) and
          H(i+1,i)*H(i,i+1) < 0. If INFO = 0 and JOB = 'E', the
          contents of H are unspecified on exit.  (The output value of
          H when INFO > 0 is given under the description of INFO
          below.) \n
 \n
          Unlike earlier versions of SHSEQR, this subroutine may
          explicitly H(i,j) = 0 for i > j and j = 1, 2, ... ILO-1
          or j = IHI+1, IHI+2, ... N. \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the array H. LDH >= fla_max(1,N). \n
 * @param[out] WR
          WR is REAL array, dimension (N) \n
 * @param[out] WI
          WI is REAL array, dimension (N) \n
          The real and imaginary parts, respectively, of the computed
          eigenvalues. If two eigenvalues are computed as a complex
          conjugate pair, they are stored in consecutive elements of
          WR and WI, say the i-th and (i+1)th, with WI(i) > 0 and
          WI(i+1) < 0. If JOB = 'S', the eigenvalues are stored in
          the same order as on the diagonal of the Schur form returned
          in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
          diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
          WI(i+1) = -WI(i). \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,N) \n
          If COMPZ = 'N', Z is not referenced. \n
          If COMPZ = 'I', on entry Z need not be set and on exit, \n
          if INFO = 0, Z contains the orthogonal matrix Z of the Schur
          vectors of H.  If COMPZ = 'V', on entry Z must contain an
          N-by-N matrix Q, which is assumed to be equal to the unit
          matrix except for the submatrix Z(ILO:IHI,ILO:IHI). \n On exit,
          if INFO = 0, Z contains Q*Z.
          Normally Q is the orthogonal matrix generated by SORGHR
          after the call to SGEHRD which formed the Hessenberg matrix
          H. (The output value of Z when INFO > 0 is given under
          the description of INFO below.) \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  if COMPZ = 'I' or
          COMPZ = 'V', then LDZ >= MAX(1,N).  Otherwise, LDZ >= 1. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (LWORK) \n
          On exit, if INFO = 0, WORK(1) returns an estimate of
          the optimal value for LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N)
          is sufficient and delivers very good and sometimes
          optimal performance.  However, LWORK as large as 11*N
          may be required for optimal performance.  A workspace
          query is recommended to determine the optimal workspace
          size. \n
\n
          If LWORK = -1, then SHSEQR does a workspace query.
          In this case, SHSEQR checks the input parameters and
          estimates the optimal workspace size for the given
          values of N, ILO and IHI.  The estimate is returned
          in WORK(1).  No error message related to LWORK is
          issued by XERBLA.  Neither H nor Z are accessed. \n
 * @param[out]	INFO
          INFO is INTEGER \n
             = 0:  successful exit \n
             < 0:  if INFO = -i, the i-th argument had an illegal
                    value \n
             > 0:  if INFO = i, SHSEQR failed to compute all of
                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
                and WI contain those eigenvalues which have been
                successfully computed.  (Failures are rare.) \n
 \n
                If INFO > 0 and JOB = 'E', then on exit, the
                remaining unconverged eigenvalues are the eigen-
                values of the upper Hessenberg matrix rows and
                columns ILO through INFO of the final, output
                value of H. \n
 \n
                If INFO > 0 and JOB   = 'S', then on exit
 \n
           (*)  (initial value of H)*U  = U*(final value of H)
 \n
                where U is an orthogonal matrix.  The final
                value of H is upper Hessenberg and quasi-triangular
                in rows and columns INFO+1 through IHI.
 \n
                If INFO > 0 and COMPZ = 'V', then on exit
 \n
                  (final value of Z)  =  (initial value of Z)*U
 \n
                where U is the orthogonal matrix in (*) (regard-
                less of the value of JOB.)
 \n
                If INFO > 0 and COMPZ = 'I', then on exit
                      (final value of Z)  = U
                where U is the orthogonal matrix in (*) (regard-
                less of the value of JOB.)
 \n
                If INFO > 0 and COMPZ = 'N', then Z is not
                accessed. \n

 *  * */
    template <typename T>
    void hseqr(char *job, char *compz, integer *n, integer *ilo, integer *ihi, T *h, integer *ldh,
               T *wr, T *wi, T *z, integer *ldz, T *work, integer *lwork, integer *info)
    {
        hseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
    }
    template <typename T>
    void hseqr(char *job, char *compz, integer *n, integer *ilo, integer *ihi, T *h, integer *ldh,
               T *w, T *z, integer *ldz, T *work, integer *lwork, integer *info)
    {
        hseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
    }
    /** @}*/ // end of hseqr

    /** @defgroup hsein hsein
     * @ingroup NYS
     * @{
     */
    /*! @brief HSEIN uses inverse iteration to find specified right and/or left \n
     eigenvectors of a real upper Hessenberg matrix H

 * @details
 * \b Purpose:
    \verbatim
     HSEIN uses inverse iteration to find specified right and/or left
     eigenvectors of a real upper Hessenberg matrix H.

     The right eigenvector x and the left eigenvector y of the matrix H
     corresponding to an eigenvalue w are defined by:

                  H * x = w * x,     y**h * H = w * y**h

     where y**h denotes the conjugate transpose of the vector y.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'R': compute right eigenvectors only; \n
          = 'L': compute left eigenvectors only; \n
          = 'B': compute both right and left eigenvectors. \n
 * @param[in] EIGSRC
          EIGSRC is CHARACTER*1 \n
          Specifies the source of eigenvalues supplied in (WR,WI): \n
          = 'Q': the eigenvalues were found using SHSEQR; thus, if
                 H has zero subdiagonal elements, and so is
                 block-triangular, then the j-th eigenvalue can be
                 assumed to be an eigenvalue of the block containing
                 the j-th row/column.  This property allows SHSEIN to
                 perform inverse iteration on just one diagonal block. \n
          = 'N': no assumptions are made on the correspondence
                 between eigenvalues and diagonal blocks.  In this
                 case, SHSEIN must always perform inverse iteration
                 using the whole matrix H. \n
 * @param[in] INITV
          INITV is CHARACTER*1 \n
          = 'N': no initial vectors are supplied; \n
          = 'U': user-supplied initial vectors are stored in the arrays
                 VL and/or VR. \n
 * @param[in,out] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          Specifies the eigenvectors to be computed. To select the
          real eigenvector corresponding to a real eigenvalue WR(j),
          SELECT(j) must be set to .TRUE.. To select the complex
          eigenvector corresponding to a complex eigenvalue
          (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)),
          either SELECT(j) or SELECT(j+1) or both must be set to
          .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is
          .FALSE.. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix H.  N >= 0. \n
 * @param[in] H
          H is REAL array, dimension (LDH,N) \n
          The upper Hessenberg matrix H.
          If a NaN is detected in H, the routine will return with INFO=-6. \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the array H.  LDH >= fla_max(1,N). \n
 * @param[in,out] WR
          WR is REAL array, dimension (N) \n
 * @param[in] WI
          WI is REAL array, dimension (N) \n
          On entry, the real and imaginary parts of the eigenvalues of
          H; a complex conjugate pair of eigenvalues must be stored in
          consecutive elements of WR and WI. \n
          On exit, WR may have been altered since close eigenvalues
          are perturbed slightly in searching for independent
          eigenvectors. \n
 * @param[in,out] VL
          VL is REAL array, dimension (LDVL,MM) \n
          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must
          contain starting vectors for the inverse iteration for the
          left eigenvectors; the starting vector for each eigenvector
          must be in the same column(s) in which the eigenvector will
          be stored. \n
          On exit, if SIDE = 'L' or 'B', the left eigenvectors
          specified by SELECT will be stored consecutively in the
          columns of VL, in the same order as their eigenvalues. A
          complex eigenvector corresponding to a complex eigenvalue is
          stored in two consecutive columns, the first holding the real
          part and the second the imaginary part.
          If SIDE = 'R', VL is not referenced. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the array VL. \n
          LDVL >= fla_max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise. \n
 * @param[in,out] VR
          VR is REAL array, dimension (LDVR,MM) \n
          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must
          contain starting vectors for the inverse iteration for the
          right eigenvectors; the starting vector for each eigenvector
          must be in the same column(s) in which the eigenvector will
          be stored. \n
          On exit, if SIDE = 'R' or 'B', the right eigenvectors
          specified by SELECT will be stored consecutively in the
          columns of VR, in the same order as their eigenvalues. A
          complex eigenvector corresponding to a complex eigenvalue is
          stored in two consecutive columns, the first holding the real
          part and the second the imaginary part.
          If SIDE = 'L', VR is not referenced. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR.
          LDVR >= fla_max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise. \n
 * @param[in] MM
          MM is INTEGER \n
          The number of columns in the arrays VL and/or VR. MM >= M. \n
 * @param[out] M
          M is INTEGER \n
          The number of columns in the arrays VL and/or VR required to
          store the eigenvectors; each selected real eigenvector
          occupies one column and each selected complex eigenvector
          occupies two columns. \n
 * @param[out] IFAILL
          IFAILL is INTEGER array, dimension (MM) \n
          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left
          eigenvector in the i-th column of VL (corresponding to the
          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the
          eigenvector converged satisfactorily. If the i-th and (i+1)th
          columns of VL hold a complex eigenvector, then IFAILL(i) and
          IFAILL(i+1) are set to the same value. \n
          If SIDE = 'R', IFAILL is not referenced. \n
 * @param[out] IFAILR
          IFAILR is INTEGER array, dimension (MM) \n
          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right
          eigenvector in the i-th column of VR (corresponding to the
          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the
          eigenvector converged satisfactorily. If the i-th and (i+1)th
          columns of VR hold a complex eigenvector, then IFAILR(i) and
          IFAILR(i+1) are set to the same value. \n
          If SIDE = 'L', IFAILR is not referenced. \n
 * @param[out]	WORK
          WORK is REAL array, dimension ((N+2)*N) \n
 * @param[out]	IFAILL
          IFAILL is INTEGER array, dimension (MM) \n
          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left
          eigenvector in the i-th column of VL (corresponding to the
          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the
          eigenvector converged satisfactorily. If the i-th and (i+1)th
          columns of VL hold a complex eigenvector, then IFAILL(i) and
          IFAILL(i+1) are set to the same value. \n
          If SIDE = 'R', IFAILL is not referenced. \n
 * @param[out]	IFAILR
          IFAILR is INTEGER array, dimension (MM) \n
          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right
          eigenvector in the i-th column of VR (corresponding to the
          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the
          eigenvector converged satisfactorily. If the i-th and (i+1)th
          columns of VR hold a complex eigenvector, then IFAILR(i) and
          IFAILR(i+1) are set to the same value. \n
          If SIDE = 'L', IFAILR is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, i is the number of eigenvectors which
                failed to converge; see IFAILL and IFAILR for further
                details. \n

 *  * */
    template <typename T>
    void hsein(char *job, char *eigsrc, char *initv, logical *select, integer *n, T *h,
               integer *ldh, T *wr, T *wi, T *vl, integer *ldvl, T *vr, integer *ldvr, integer *mm,
               integer *m, T *work, integer *ifaill, integer *ifailr, integer *info)
    {
        hsein(job, eigsrc, initv, select, n, h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work,
              ifaill, ifailr, info);
    }
    template <typename T, typename Ta>
    void hsein(char *job, char *eigsrc, char *initv, logical *select, integer *n, T *h,
               integer *ldh, T *w, T *vl, integer *ldvl, T *vr, integer *ldvr, integer *mm,
               integer *m, T *work, Ta *rwork, integer *ifaill, integer *ifailr, integer *info)
    {
        hsein(job, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, rwork,
              ifaill, ifailr, info);
    }
    /** @}*/ // end of hsein

    /** @defgroup trevc trevc
     * @ingroup NYS
     * @{
     */
    /*! @brief TREVC computes some or all of the right and/or left eigenvectors of  \n
     a real upper quasi-triangular matrix T.

 * @details
 * \b Purpose:
    \verbatim
     TREVC computes some or all of the right and/or left eigenvectors of
     a real upper quasi-triangular matrix T.
     Matrices of this type are produced by the Schur factorization of
     a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR.

     The right eigenvector x and the left eigenvector y of T corresponding
     to an eigenvalue w are defined by:

        T*x = w*x,     (y**H)*T = w*(y**H)

     where y**H denotes the conjugate transpose of y.
     The eigenvalues are not input to this routine, but are read directly
     from the diagonal blocks of T.

     This routine returns the matrices X and/or Y of right and left
     eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     input matrix.  If Q is the orthogonal factor that reduces a matrix
     A to Schur form T, then Q*X and Q*Y are the matrices of right and
     left eigenvectors of A.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'R':  compute right eigenvectors only; \n
          = 'L':  compute left eigenvectors only; \n
          = 'B':  compute both right and left eigenvectors. \n
 * @param[in] HOWMNY
          HOWMNY is CHARACTER*1 \n
          = 'A':  compute all right and/or left eigenvectors; \n
          = 'B':  compute all right and/or left eigenvectors,
                  backtransformed by the matrices in VR and/or VL; \n
          = 'S':  compute selected right and/or left eigenvectors,
                  as indicated by the logical array SELECT. \n
 * @param[in,out] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
          computed. \n
          If w(j) is a real eigenvalue, the corresponding real
          eigenvector is computed if SELECT(j) is .TRUE.. \n
          If w(j) and w(j+1) are the real and imaginary parts of a
          complex eigenvalue, the corresponding complex eigenvector is
          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and
          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to
          .FALSE.. \n
          Not referenced if HOWMNY = 'A' or 'B'. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. N >= 0. \n
 * @param[in] T
          T is REAL array, dimension (LDT,N) \n
          The upper quasi-triangular matrix T in Schur canonical form. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= fla_max(1,N). \n
 * @param[in,out] VL
          VL is REAL array, dimension (LDVL,MM) \n
          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
          contain an N-by-N matrix Q (usually the orthogonal matrix Q
          of Schur vectors returned by SHSEQR). \n
          On exit, if SIDE = 'L' or 'B', VL contains: \n
          if HOWMNY = 'A', the matrix Y of left eigenvectors of T; \n
          if HOWMNY = 'B', the matrix Q*Y; \n
          if HOWMNY = 'S', the left eigenvectors of T specified by
                           SELECT, stored consecutively in the columns
                           of VL, in the same order as their
                           eigenvalues. \n
          A complex eigenvector corresponding to a complex eigenvalue
          is stored in two consecutive columns, the first holding the
          real part, and the second the imaginary part.
          Not referenced if SIDE = 'R'. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the array VL.  LDVL >= 1, and if
          SIDE = 'L' or 'B', LDVL >= N. \n
 * @param[in,out] VR
          VR is REAL array, dimension (LDVR,MM) \n
          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
          contain an N-by-N matrix Q (usually the orthogonal matrix Q
          of Schur vectors returned by SHSEQR). \n
          On exit, if SIDE = 'R' or 'B', VR contains: \n
          if HOWMNY = 'A', the matrix X of right eigenvectors of T; \n
          if HOWMNY = 'B', the matrix Q*X; \n
          if HOWMNY = 'S', the right eigenvectors of T specified by
                           SELECT, stored consecutively in the columns
                           of VR, in the same order as their
                           eigenvalues. \n
          A complex eigenvector corresponding to a complex eigenvalue
          is stored in two consecutive columns, the first holding the
          real part and the second the imaginary part. \n
          Not referenced if SIDE = 'L'. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR.  LDVR >= 1, and if
          SIDE = 'R' or 'B', LDVR >= N. \n
 * @param[in] MM
          MM is INTEGER \n
          The number of columns in the arrays VL and/or VR. MM >= M. \n
 * @param[out] M
          M is INTEGER \n
          The number of columns in the arrays VL and/or VR actually
          used to store the eigenvectors. \n
          If HOWMNY = 'A' or 'B', M is set to N.
          Each selected real eigenvector occupies one column and each
          selected complex eigenvector occupies two columns. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void trevc(char *side, char *howmny, logical *select, integer *n, T *t, integer *ldt, T *vl,
               integer *ldvl, T *vr, integer *ldvr, integer *mm, integer *m, T *work, integer *info)
    {
        trevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info);
    }
    template <typename T, typename Ta>
    void trevc(char *side, char *howmny, logical *select, integer *n, T *t, integer *ldt, T *vl,
               integer *ldvl, T *vr, integer *ldvr, integer *mm, integer *m, T *work, Ta *rwork,
               integer *info)
    {
        trevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
    }
    /** @}*/ // end of trevc

    /** @defgroup trevc3 trevc3
     * @ingroup NYS
     * @{
     */
    /*! @brief TREVC3 computes some or all of the right and/or left eigenvectors \n
     of a real upper quasi-triangular matrix T

 * @details
 * \b Purpose:
    \verbatim
    TREVC3 computes some or all of the right and/or left eigenvectors of
    a complex upper triangular matrix T.
    Matrices of this type are produced by the Schur factorization of
    a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR.

    The right eigenvector x and the left eigenvector y of T corresponding
    to an eigenvalue w are defined by:

                 T*x = w*x,    (y**H)*T = w*(y**H)

    where y**H denotes the conjugate transpose of the vector y.
    The eigenvalues are not input to this routine, but are read directly
    from the diagonal of T.

    This routine   returns the matrices X and/or Y of right and left
    eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
    input matrix. If Q is the unitary factor that reduces a matrix A to
    Schur form T, then Q*X and Q*Y are the matrices of right and left
    eigenvectors of A.

    This uses a Level 3 BLAS version of the back transformation.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'R':  compute right eigenvectors only; \n
          = 'L':  compute left eigenvectors only; \n
          = 'B':  compute both right and left eigenvectors. \n
 * @param[in] HOWMNY
          HOWMNY is CHARACTER*1 \n
          = 'A':  compute all right and/or left eigenvectors; \n
          = 'B':  compute all right and/or left eigenvectors,
                  backtransformed using the matrices supplied in
                  VR and/or VL; \n
          = 'S':  compute selected right and/or left eigenvectors,
                  as indicated by the logical array SELECT. \n
 * @param[in] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
          computed. \n
          The eigenvector corresponding to the j-th eigenvalue is
          computed if SELECT(j) = .TRUE..
          Not referenced if HOWMNY = 'A' or 'B'. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. N >= 0. \n
 * @param[in,out] T
          T is COMPLEX array, dimension (LDT,N) \n
          The upper triangular matrix T.  T is modified, but restored
          on exit. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= fla_max(1,N). \n
 * @param[in,out] VL
          VL is COMPLEX array, dimension (LDVL,MM) \n
          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
          contain an N-by-N matrix Q (usually the unitary matrix Q of
          Schur vectors   returned by CHSEQR). \n
          On exit, if SIDE = 'L' or 'B', VL contains: \n
          if HOWMNY = 'A', the matrix Y of left eigenvectors of T; \n
          if HOWMNY = 'B', the matrix Q*Y; \n
          if HOWMNY = 'S', the left eigenvectors of T specified by
                           SELECT, stored consecutively in the columns
                           of VL, in the same order as their
                           eigenvalues. \n
          Not referenced if SIDE = 'R'. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the array VL. \n
          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N. \n
 * @param[in,out] VR
          VR is COMPLEX array, dimension (LDVR,MM) \n
          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
          contain an N-by-N matrix Q (usually the unitary matrix Q of
          Schur vectors   returned by CHSEQR). \n
          On exit, if SIDE = 'R' or 'B', VR contains: \n
          if HOWMNY = 'A', the matrix X of right eigenvectors of T; \n
          if HOWMNY = 'B', the matrix Q*X; \n
          if HOWMNY = 'S', the right eigenvectors of T specified by
                           SELECT, stored consecutively in the columns
                           of VR, in the same order as their
                           eigenvalues. \n
          Not referenced if SIDE = 'L'. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR. \n
          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N. \n
 * @param[in] MM
          MM is INTEGER \n
          The number of columns in the arrays VL and/or VR. MM >= M. \n
 * @param[out] M
          M is INTEGER \n
          The number of columns in the arrays VL and/or VR actually
          used to store the eigenvectors. \n
          If HOWMNY = 'A' or 'B', M is set to N.
          Each selected eigenvector occupies one column. \n
 * @param[out] WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of array WORK. LWORK >= fla_max(1,2*N). \n
          For optimum performance, LWORK >= N + 2*N*NB, where NB is
          the optimal blocksize. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] RWORK
          RWORK is REAL array, dimension (LRWORK) \n
 * @param[in] LRWORK
          LRWORK is INTEGER \n
          The dimension of array RWORK. LRWORK >= fla_max(1,N). \n
 \n
          If LRWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the RWORK array,   returns
          this value as the first entry of the RWORK array, and no error
          message related to LRWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void trevc3(char *side, char *howmny, logical *select, integer *n, T *t, integer *ldt, T *vl,
                integer *ldvl, T *vr, integer *ldvr, integer *mm, integer *m, T *work,
                integer *lwork, integer *info)
    {
        trevc3(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info);
    }
    template <typename T, typename Ta>
    void trevc3(char *side, char *howmny, logical *select, integer *n, T *t, integer *ldt, T *vl,
                integer *ldvl, T *vr, integer *ldvr, integer *mm, integer *m, T *work,
                integer *lwork, Ta *rwork, integer *lrwork, integer *info)
    {
        trevc3(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, rwork,
               lrwork, info);
    }
    /** @}*/ // end of trevc3

    /** @defgroup laln2 laln2
     * @ingroup NYS
     * @{
     */
    /*! @brief LALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form

 * @details
 * \b Purpose:
    \verbatim
    LALN2 solves a system of the form  (ca A - w D) X = s B
    or (ca A**T - w D) X = s B   with possible scaling ("s") and
    perturbation of A.  (A**T means A-transpose.)

    A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
    real diagonal matrix, w is a real or complex value, and X and B are
    NA x 1 matrices -- real if w is real, complex if w is complex.  NA
    may be 1 or 2.

    If w is complex, X and B are represented as NA x 2 matrices,
    the first column of each being the real part and the second
    being the imaginary part.

    "s" is a scaling factor (<= 1), computed by SLALN2, which is
    so chosen that X can be computed without overflow.  X is further
    scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
    than overflow.

    If both singular values of (ca A - w D) are less than SMIN,
    SMIN*identity will be used instead of (ca A - w D).  If only one
    singular value is less than SMIN, one element of (ca A - w D) will be
    perturbed enough to make the smallest singular value roughly SMIN.
    If both singular values are at least SMIN, (ca A - w D) will not be
    perturbed.  In any case, the perturbation will be at most some small
    multiple of fla_max(SMIN, ulp*norm(ca A - w D)).  The singular values
    are computed by infinity-norm approximations, and thus will only be
    correct to a factor of 2 or so.

    Note: all input quantities are assumed to be smaller than overflow
    by a reasonable factor.  (See BIGNUM.)
    \endverbatim

 * @param[in] LTRANS
          LTRANS is LOGICAL \n
          =.TRUE.:  A-transpose will be used. \n
          =.FALSE.: A will be used (not transposed.) \n
 * @param[in] NA
          NA is INTEGER \n
          The size of the matrix A.  It may (only) be 1 or 2. \n
 * @param[in] NW
          NW is INTEGER \n
          1 if "w" is real, 2 if "w" is complex.  It may only be 1
          or 2. \n
 * @param[in] SMIN
          SMIN is REAL \n
          The desired lower bound on the singular values of A.  This
          should be a safe distance away from underflow or overflow,
          say, between (underflow/machine precision) and  (machine
          precision * overflow).  (See BIGNUM and ULP.) \n
 * @param[in] CA
          CA is REAL \n
          The coefficient c, which A is multiplied by. \n
 * @param[in] A
          A is REAL array, dimension (LDA,NA) \n
          The NA x NA matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A.  It must be at least NA. \n
 * @param[in] D1
          D1 is REAL \n
          The 1,1 element in the diagonal matrix D. \n
 * @param[in] D2
          D2 is REAL \n
          The 2,2 element in the diagonal matrix D.  Not used if NA=1. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NW) \n
          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
          complex), column 1 contains the real part of B and column 2
          contains the imaginary part. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B.  It must be at least NA. \n
 * @param[in] WR
          WR is REAL \n
          The real part of the scalar "w". \n
 * @param[in] WI
          WI is REAL \n
          The imaginary part of the scalar "w".  Not used if NW=1. \n
 * @param[out] X
          X is REAL array, dimension (LDX,NW) \n
          The NA x NW matrix X (unknowns), as computed by SLALN2.
          If NW=2 ("w" is complex), on exit, column 1 will contain
          the real part of X and column 2 will contain the imaginary
          part. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of X.  It must be at least NA. \n
 * @param[out] SCALE
          SCALE is REAL \n
          The scale factor that B must be multiplied by to insure
          that overflow does not occur when computing X.  Thus,
          (ca A - w D) X  will be SCALE*B, not B (ignoring
          perturbations of A.)  It will be at most 1. \n
 * @param[out] XNORM
          XNORM is REAL \n
          The infinity-norm of X, when X is regarded as an NA x NW
          real matrix. \n
 * @param[out] INFO
          INFO is INTEGER \n
          An error flag.  It will be set to zero if no error occurs,
          a negative number if an argument is in error, or a positive
          number if  ca A - w D  had to be perturbed.
          The possible values are: \n
          = 0: No error occurred, and (ca A - w D) did not have to be
                 perturbed. \n
          = 1: (ca A - w D) had to be perturbed to make its smallest
               (or only) singular value greater than SMIN. \n
          NOTE: In the interests of speed, this routine does not
                check the inputs for errors. \n

 * */
    template <typename T>
    void laln2(logical *ltrans, integer *na, integer *nw, float *smin, float *ca, float *a,
               integer *lda, float *d1, float *d2, float *b, integer *ldb, float *wr, float *wi,
               float *x, integer *ldx, float *scale, float *xnorm, integer *info)
    {
        laln2(ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale, xnorm, info);
    }
    /** @}*/ // end of laln2

    /** @defgroup trsyl trsyl
     * @ingroup NYS
     * @{
     */
    /*! @brief Solving Sylvester matrix equation
    *
    * @details
    * \b Purpose:
    * \verbatim
        Solution for real Sylvester matrix equation:
            op(A)*X + X*op(B) = scale*C or
            op(A)*X - X*op(B) = scale*C,

        where op(A) = A or A**T, and  a and b are both upper quasi- triangular.
        A is M-by-M and B is n-by-n; the right hand side C and the solution X are m-by-n;
        and scale is an output scale factor, set <= 1 to avoid overflow in X.

        a and b must be in Schur canonical form (as returned by SHSEQR), that is, block upper
        triangular with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has its
        diagonal elements equal and its off-diagonal elements of opposite sign.
    \endverbatim

    * @param[in] transa
              transa is char* \n
              Specifies the option op(A): \n
              = 'N': op(A) = A (No transpose) \n
              = 'T': op(A) = A**T (Transpose) \n
              = 'C': op(A) = A**H (Conjugate transpose = Transpose) \n
    * @param[in] transb
              transb is char* \n
              Specifies the option op(B): \n
              = 'N': op(B) = B (No transpose) \n
              = 'T': op(B) = B**T (Transpose) \n
              = 'C': op(B) = B**H (Conjugate transpose = Transpose) \n
    * @param[in] isgn
              isgn is integer* \n
              Specifies the sign in the equation: \n
              = +1: solve op(A)*X + X*op(B) = scale*C \n
              = -1: solve op(A)*X - X*op(B) = scale*C \n
    * @param[in] m
              m is integer* \n
              The order of the matrix a, and the number of rows in the
              matrices X and C. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix b, and the number of columns in the
              matrices X and C. n >= 0. \n
    * @param[in] a
              a is float/double array, dimension (lda,m) \n
              The upper quasi-triangular matrix a, in Schur canonical form. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. lda >= fla_max(1,m). \n
    * @param[in] b
              b is float/double array, dimension (ldb,n) \n
              The upper quasi-triangular matrix b, in Schur canonical form. \n
    * @param[in] ldb
              ldb is integer* \n
              The leading dimension of the array b. ldb >= fla_max(1,n). \n
    * @param[in,out] c
              c is float/double array, dimension (ldc,n) \n
              On entry, the m-by-n right hand side matrix c. \n
              On exit, c is overwritten by the solution matrix X. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. ldc >= fla_max(1,m) \n
    * @param[out] scale
              scale is float/double* \n
              The scale factor, scale, set <= 1 to avoid overflow in X \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -i, the i-th argument had an illegal value \n
              = 1: A and B have common or very close eigenvalues; perturbed
                   values were used to solve the equation (but the matrices
                   A and B are unchanged). \n

    *     *  */
    template <typename T>
    void trsyl(char *transa, char *transb, integer *isgn, integer *m, integer *n, T *a,
               integer *lda, T *b, integer *ldb, T *c, integer *ldc, T *scale, integer *info)
    {
        trsyl(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
    }
    template <typename T, typename Ta>
    void trsyl(char *transa, char *transb, integer *isgn, integer *m, integer *n, T *a,
               integer *lda, T *b, integer *ldb, T *c, integer *ldc, Ta *scale, integer *info)
    {
        trsyl(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
    }
    /** @}*/ // end of trsyl

    /** @defgroup trsna trsna
     * @ingroup NYS
     * @{
     */
    /*! @brief TRSNA estimates reciprocal condition numbers for specified eigenvalues

 * @details
 * \b Purpose:
    \verbatim
     TRSNA estimates reciprocal condition numbers for specified
     eigenvalues and/or right eigenvectors of a real upper
     quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
     orthogonal).

     T must be in Schur canonical form (as returned by SHSEQR), that is,
     block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     2-by-2 diagonal block has its diagonal elements equal and its
     off-diagonal elements of opposite sign.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          Specifies whether condition numbers are required for
          eigenvalues (S) or eigenvectors (SEP): \n
          = 'E': for eigenvalues only (S); \n
          = 'V': for eigenvectors only (SEP); \n
          = 'B': for both eigenvalues and eigenvectors (S and SEP). \n
 * @param[in] HOWMNY
          HOWMNY is CHARACTER*1 \n
          = 'A': compute condition numbers for all eigenpairs; \n
          = 'S': compute condition numbers for selected eigenpairs
                 specified by the array SELECT. \n
 * @param[in] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
          condition numbers are required. To select condition numbers
          for the eigenpair corresponding to a real eigenvalue w(j),
          SELECT(j) must be set to .TRUE.. To select condition numbers
          corresponding to a complex conjugate pair of eigenvalues w(j)
          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
          set to .TRUE.. \n
          If HOWMNY = 'A', SELECT is not referenced. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. N >= 0. \n
 * @param[in] T
          T is REAL array, dimension (LDT,N) \n
          The upper quasi-triangular matrix T, in Schur canonical form. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= fla_max(1,N). \n
 * @param[in] VL
          VL is REAL array, dimension (LDVL,M) \n
          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
          must be stored in consecutive columns of VL, as returned by
          SHSEIN or STREVC. \n
          If JOB = 'V', VL is not referenced. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the array VL.
          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N. \n
 * @param[in] VR
          VR is REAL array, dimension (LDVR,M) \n
          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
          must be stored in consecutive columns of VR, as returned by
          SHSEIN or STREVC. \n
          If JOB = 'V', VR is not referenced. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR. \n
          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N. \n
 * @param[out] S
          S is REAL array, dimension (MM) \n
          If JOB = 'E' or 'B', the reciprocal condition numbers of the
          selected eigenvalues, stored in consecutive elements of the
          array. For a complex conjugate pair of eigenvalues two
          consecutive elements of S are set to the same value. Thus
          S(j), SEP(j), and the j-th columns of VL and VR all
          correspond to the same eigenpair (but not in general the
          j-th eigenpair, unless all eigenpairs are selected). \n
          If JOB = 'V', S is not referenced. \n
 * @param[out] SEP
          SEP is REAL array, dimension (MM) \n
          If JOB = 'V' or 'B', the estimated reciprocal condition
          numbers of the selected eigenvectors, stored in consecutive
          elements of the array. For a complex eigenvector two
          consecutive elements of SEP are set to the same value. If
          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)
          is set to 0; this can only occur when the true value would be
          very small anyway. \n
          If JOB = 'E', SEP is not referenced. \n
 * @param[in] MM
          MM is INTEGER \n
          The number of elements in the arrays S (if JOB = 'E' or 'B')
          and/or SEP (if JOB = 'V' or 'B'). MM >= M. \n
 * @param[out] M
          M is INTEGER \n
          The number of elements of the arrays S and/or SEP actually
          used to store the estimated condition numbers.
          If HOWMNY = 'A', M is set to N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (LDWORK,N+6) \n
          If JOB = 'E', WORK is not referenced. \n
 * @param[in]	LDWORK
          LDWORK is INTEGER \n
          The leading dimension of the array WORK.
          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (2*(N-1)) \n
          If JOB = 'E', IWORK is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void trsna(char *job, char *howmny, logical *select, integer *n, T *t, integer *ldt, T *vl,
               integer *ldvl, T *vr, integer *ldvr, T *s, T *sep, integer *mm, integer *m, T *work,
               integer *ldwork, integer *iwork, integer *info)
    {
        trsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, m, work, ldwork,
              iwork, info);
    }
    template <typename T, typename Ta>
    void trsna(char *job, char *howmny, logical *select, integer *n, T *t, integer *ldt, T *vl,
               integer *ldvl, T *vr, integer *ldvr, Ta *s, Ta *sep, integer *mm, integer *m,
               T *work, integer *ldwork, Ta *rwork, integer *info)
    {
        trsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, m, work, ldwork,
              rwork, info);
    }
    /** @}*/ // end of trsna

    /** @defgroup laqtr laqtr
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQTR solves a real quasi-triangular system of equations, or a complex \n
     quasi-triangular system of special form, in real arithmetic
 * @details
 * \b Purpose:
    \verbatim
    LAQTR solves the real quasi-triangular system

                 op(T)*p = scale*c,              if LREAL = .TRUE.

    or the complex quasi-triangular systems

               op(T + iB)*(p+iq) = scale*(c+id), if LREAL = .FALSE.

    in real arithmetic, where T is upper quasi-triangular.
    If LREAL = .FALSE., then the first diagonal block of T must be
    1 by 1, B is the specially structured matrix

                   B = [ b(1) b(2) ... b(n) ]
                       [       w            ]
                       [           w        ]
                       [              .     ]
                       [                 w  ]

    op(A) = A or A**T, A**T denotes the transpose of
    matrix A.

    On input, X = [ c ].  On output, X = [ p ].
                  [ d ]                  [ q ]

    This subroutine is designed for the condition number estimation
    in routine STRSNA.
    \endverbatim

 * @param[in] LTRAN
          LTRAN is LOGICAL \n
          On entry, LTRAN specifies the option of conjugate transpose: \n
             = .FALSE.,   op(T+i*B) = T+i*B, \n
             = .TRUE.,    op(T+i*B) = (T+i*B)**T. \n
 * @param[in] LREAL
          LREAL is LOGICAL \n
          On entry, LREAL specifies the input matrix structure: \n
             = .FALSE.,   the input is complex \n
             = .TRUE.,    the input is real \n
 * @param[in] N
          N is INTEGER \n
          On entry, N specifies the order of T+i*B. N >= 0. \n
 * @param[in] T
          T is REAL array, dimension (LDT,N) \n
          On entry, T contains a matrix in Schur canonical form. \n
          If LREAL = .FALSE., then the first diagonal block of T must
          be 1 by 1. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the matrix T. LDT >= fla_max(1,N). \n
 * @param[in] B
          B is REAL array, dimension (N) \n
          On entry, B contains the elements to form the matrix
          B as described above. \n
          If LREAL = .TRUE., B is not referenced. \n
 * @param[in] W
          W is REAL \n
          On entry, W is the diagonal element of the matrix B.
          If LREAL = .TRUE., W is not referenced. \n
 * @param[out] SCALE
          SCALE is REAL \n
          On exit, SCALE is the scale factor. \n
 * @param[in,out] X
          X is REAL array, dimension (2*N) \n
          On entry, X contains the right hand side of the system.
          On exit, X is overwritten by the solution. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          On exit, INFO is set to \n
             0: successful exit. \n
               1: the some diagonal 1 by 1 block has been perturbed by
                  a small number SMIN to keep nonsingularity. \n
               2: the some diagonal 2 by 2 block has been perturbed by
                  a small number in SLALN2 to keep nonsingularity. \n
          NOTE: In the interests of speed, this routine does not
                check the inputs for errors.  \n

 *  * */
    template <typename T>
    void laqtr(logical *ltran, logical *lreal, integer *n, T *t, integer *ldt, T *b, T *w, T *scale,
               T *x, T *work, integer *info)
    {
        laqtr(ltran, lreal, n, t, ldt, b, w, scale, x, work, info);
    }
    /** @}*/ // end of laqtr

    /** @defgroup trexc trexc
     * @ingroup NYS
     * @{
     */
    /*! @brief TREXC reorders the real Schur factorization of a real matrix A = Q*T*Q**T

 * @details
 * \b Purpose:
    \verbatim
     TREXC reorders the real Schur factorization of a real matrix
     A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
     moved to row ILST.

     The real Schur form T is reordered by an orthogonal similarity
     transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
     is updated by postmultiplying it with Z.

     T must be in Schur canonical form (as returned by SHSEQR), that is,
     block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     2-by-2 diagonal block has its diagonal elements equal and its
     off-diagonal elements of opposite sign.
    \endverbatim

 * @param[in] COMPQ
          COMPQ is CHARACTER*1 \n
          = 'V':  update the matrix Q of Schur vectors; \n
          = 'N':  do not update Q. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. N >= 0.
          If N == 0 arguments ILST and IFST may be any value. \n
 * @param[in,out] T
          T is REAL array, dimension (LDT,N) \n
          On entry, the upper quasi-triangular matrix T, in Schur
          Schur canonical form. \n
          On exit, the reordered upper quasi-triangular matrix, again
          in Schur canonical form. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= fla_max(1,N). \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
          On exit, if COMPQ = 'V', Q has been postmultiplied by the
          orthogonal transformation matrix Z which reorders T.
          If COMPQ = 'N', Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= 1, and if
          COMPQ = 'V', LDQ >= fla_max(1,N). \n
 * @param[in,out] IFST
          IFST is INTEGER \n
 * @param[in,out] ILST
          ILST is INTEGER \n
          Specify the reordering of the diagonal blocks of T.
          The block with row index IFST is moved to row ILST, by a
          sequence of transpositions between adjacent blocks.
          On exit, if IFST pointed on entry to the second row of a
          2-by-2 block, it is changed to point to the first row; ILST
          always points to the first row of the block in its final
          position (which may differ from its input value by +1 or -1).
          1 <= IFST <= N; 1 <= ILST <= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          = 1:  two adjacent blocks were too close to swap (the problem
                is very ill-conditioned); T may have been partially
                reordered, and ILST points to the first row of the
                current position of the block being moved. \n

 *  * */
    template <typename T>
    void trexc(char *compq, integer *n, T *t, integer *ldt, T *q, integer *ldq, integer *ifst,
               integer *ilst, T *work, integer *info)
    {
        trexc(compq, n, t, ldt, q, ldq, ifst, ilst, work, info);
    }
    template <typename T>
    void trexc(char *compq, integer *n, T *t, integer *ldt, T *q, integer *ldq, integer *ifst,
               integer *ilst, integer *info)
    {
        trexc(compq, n, t, ldt, q, ldq, ifst, ilst, info);
    }
    /** @}*/ // end of trexc

    /** @defgroup trsen trsen
     * @ingroup NYS
     * @{
     */
    /*! @brief TRSEN reorders the real Schur factorization of a real matrix A = Q*T*Q**T

 * @details
 * \b Purpose:
    \verbatim
     TRSEN reorders the real Schur factorization of a real matrix
     A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
     the leading diagonal blocks of the upper quasi-triangular matrix T,
     and the leading columns of Q form an orthonormal basis of the
     corresponding right invariant subspace.

     Optionally the routine computes the reciprocal condition numbers of
     the cluster of eigenvalues and/or the invariant subspace.

     T must be in Schur canonical form (as returned by SHSEQR), that is,
     block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     2-by-2 diagonal block has its diagonal elements equal and its
     off-diagonal elements of opposite sign.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          Specifies whether condition numbers are required for the
          cluster of eigenvalues (S) or the invariant subspace (SEP): \n
          = 'N': none; \n
          = 'E': for eigenvalues only (S); \n
          = 'V': for invariant subspace only (SEP); \n
          = 'B': for both eigenvalues and invariant subspace (S and
                 SEP). \n
 * @param[in] COMPQ
          COMPQ is CHARACTER*1 \n
          = 'V': update the matrix Q of Schur vectors; \n
          = 'N': do not update Q. \n
 * @param[in] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          SELECT specifies the eigenvalues in the selected cluster. To
          select a real eigenvalue w(j), SELECT(j) must be set to
          .TRUE.. To select a complex conjugate pair of eigenvalues
          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
          either SELECT(j) or SELECT(j+1) or both must be set to
          .TRUE.; a complex conjugate pair of eigenvalues must be
          either both included in the cluster or both excluded. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. N >= 0. \n
 * @param[in,out] T
          T is REAL array, dimension (LDT,N) \n
          On entry, the upper quasi-triangular matrix T, in Schur
          canonical form. \n
          On exit, T is overwritten by the reordered matrix T, again in
          Schur canonical form, with the selected eigenvalues in the
          leading diagonal blocks. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= fla_max(1,N). \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. \n
          On exit, if COMPQ = 'V', Q has been postmultiplied by the
          orthogonal transformation matrix which reorders T; the
          leading M columns of Q form an orthonormal basis for the
          specified invariant subspace. \n
          If COMPQ = 'N', Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.
          LDQ >= 1; and if COMPQ = 'V', LDQ >= N. \n
 * @param[out] WR
          WR is REAL array, dimension (N) \n
 * @param[out] WI
          WI is REAL array, dimension (N) \n
          The real and imaginary parts, respectively, of the reordered
          eigenvalues of T. The eigenvalues are stored in the same
          order as on the diagonal of T, with WR(i) = T(i,i) and, if
          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
          sufficiently ill-conditioned, then its value may differ
          significantly from its value before reordering. \n
 * @param[out] M
          M is INTEGER \n
          The dimension of the specified invariant subspace. \n
          0 < = M <= N. \n
 * @param[out] S
          S is REAL \n
          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
          condition number for the selected cluster of eigenvalues.
          S cannot underestimate the true reciprocal condition number
          by more than a factor of sqrt(N). If M = 0 or N, S = 1. \n
          If JOB = 'N' or 'V', S is not referenced. \n
 * @param[out] SEP
          SEP is REAL \n
          If JOB = 'V' or 'B', SEP is the estimated reciprocal
          condition number of the specified invariant subspace. If
          M = 0 or N, SEP = norm(T). \n
          If JOB = 'N' or 'E', SEP is not referenced. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If JOB = 'N', LWORK >= fla_max(1,N); \n
          if JOB = 'E', LWORK >= fla_max(1,M*(N-M)); \n
          if JOB = 'V' or 'B', LWORK >= fla_max(1,2*M*(N-M)). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          If JOB = 'N' or 'E', LIWORK >= 1; \n
          if JOB = 'V' or 'B', LIWORK >= fla_max(1,M*(N-M)). \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the IWORK array,
          returns this value as the first entry of the IWORK array, and
          no error message related to LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          = 1: reordering of T failed because some eigenvalues are too
               close to separate (the problem is very ill-conditioned);
               T may have been partially reordered, and WR and WI
               contain the eigenvalues in the same order as in T; S and
               SEP (if requested) are set to zero. \n

 *  * */
    template <typename T>
    void trsen(char *job, char *compq, logical *select, integer *n, T *t, integer *ldt, T *q,
               integer *ldq, T *wr, T *wi, integer *m, T *s, T *sep, T *work, integer *lwork,
               integer *iwork, integer *liwork, integer *info)
    {
        trsen(job, compq, select, n, t, ldt, q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork,
              info);
    }
    template <typename T, typename Ta>
    void trsen(char *job, char *compq, logical *select, integer *n, T *t, integer *ldt, T *q,
               integer *ldq, T *w, integer *m, Ta *s, Ta *sep, T *work, integer *lwork,
               integer *info)
    {
        trsen(job, compq, select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info);
    }
    /** @}*/ // end of trsen

    /** @defgroup laexc laexc
     * @ingroup NYS
     * @{
     */
    /*! @brief LAEXC swaps adjacent diagonal blocks of a real upper quasi-triangular \n
     matrix in Schur canonical form, by an orthogonal similarity transformation
 * @details
 * \b Purpose:
    \verbatim
    LAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
    an upper quasi-triangular matrix T by an orthogonal similarity
    transformation.

    T must be in Schur canonical form, that is, block upper triangular
    with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
    has its diagonal elemnts equal and its off-diagonal elements of
    opposite sign.
    \endverbatim

 * @param[in] WANTQ
          WANTQ is LOGICAL \n
          = .TRUE. : accumulate the transformation in the matrix Q; \n
          = .FALSE.: do not accumulate the transformation. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. N >= 0. \n
 * @param[in,out] T
          T is REAL array, dimension (LDT,N) \n
          On entry, the upper quasi-triangular matrix T, in Schur
          canonical form. \n
          On exit, the updated matrix T, again in Schur canonical form. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T. LDT >= fla_max(1,N). \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
          On exit, if WANTQ is .TRUE., the updated matrix Q. \n
          If WANTQ is .FALSE., Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.
          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N. \n
 * @param[in] J1
          J1 is INTEGER \n
          The index of the first row of the first block T11. \n
 * @param[in] N1
          N1 is INTEGER \n
          The order of the first block T11. N1 = 0, 1 or 2. \n
 * @param[in] N2
          N2 is INTEGER \n
          The order of the second block T22. N2 = 0, 1 or 2. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          = 1: the transformed matrix T would be too far from Schur
               form; the blocks are not swapped and T and Q are
               unchanged. \n

 * */
    template <typename T>
    void laexc(logical *wantq, integer *n, T *t, integer *ldt, T *q, integer *ldq, integer *j1,
               integer *n1, integer *n2, T *work, integer *info)
    {
        laexc(wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info);
    }
    /** @}*/ // end of laexc

    /** @defgroup lanv2 lanv2
     * @ingroup NYS
     * @{
     */
    /*! @brief LANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in
 standard form

 * @details
 * \b Purpose:
    \verbatim
    LANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
    matrix in standard form:

         [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
         [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]

    where either
    1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
    2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
    conjugate eigenvalues.
    \endverbatim

 * @param[in,out] A
          A is REAL \n
 * @param[in,out] B
          B is REAL \n
 * @param[in,out] C
          C is REAL \n
 * @param[in,out] D
          D is REAL \n
          On entry, the elements of the input matrix.
          On exit, they are overwritten by the elements of the
          standardised Schur form. \n
 * @param[out] RT1R
          RT1R is REAL \n
 * @param[out] RT1I
          RT1I is REAL \n
 * @param[out] RT2R
          RT2R is REAL \n
 * @param[out] RT2I
          RT2I is REAL \n
          The real and imaginary parts of the eigenvalues. If the
          eigenvalues are a complex conjugate pair, RT1I > 0. \n
 * @param[out] CS
          CS is REAL \n
 * @param[out] SN
          SN is REAL \n
          Parameters of the rotation matrix. \n

 * */
    template <typename T>
    void lanv2(T *a, T *b, T *c, T *d, T *rt1r, T *rt1i, T *rt2r, T *rt2i, T *cs, T *sn)
    {
        lanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn);
    }
    /** @}*/ // end of lanv2

    /** @defgroup laein laein
     * @ingroup NYS
     * @{
     */
    /*! @brief LAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by
 inverse iteration

 * @details
 * \b Purpose:
    \verbatim
    LAEIN uses inverse iteration to find a right or left eigenvector
    corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
    matrix H.
    \endverbatim

 * @param[in] RIGHTV
          RIGHTV is LOGICAL \n
          = .TRUE. : compute right eigenvector; \n
          = .FALSE.: compute left eigenvector. \n
 * @param[in] NOINIT
          NOINIT is LOGICAL \n
          = .TRUE. : no initial vector supplied in (VR,VI). \n
          = .FALSE.: initial vector supplied in (VR,VI). \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix H.  N >= 0. \n
 * @param[in] H
          H is REAL array, dimension (LDH,N) \n
          The upper Hessenberg matrix H. \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the array H.  LDH >= fla_max(1,N). \n
 * @param[in] WR
          WR is REAL \n
 * @param[in] WI
          WI is REAL \n
          The real and imaginary parts of the eigenvalue of H whose
          corresponding right or left eigenvector is to be computed. \n
 * @param[in,out] VR
          VR is REAL array, dimension (N) \n
 * @param[in,out] VI
          VI is REAL array, dimension (N) \n
          On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain
          a real starting vector for inverse iteration using the real
          eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI
          must contain the real and imaginary parts of a complex
          starting vector for inverse iteration using the complex
          eigenvalue (WR,WI); otherwise VR and VI need not be set.
          On exit, if WI = 0.0 (real eigenvalue), VR contains the
          computed real eigenvector; if WI.ne.0.0 (complex eigenvalue),
          VR and VI contain the real and imaginary parts of the
          computed complex eigenvector. The eigenvector is normalized
          so that the component of largest magnitude has magnitude 1;
          here the magnitude of a complex number (x,y) is taken to be
          |x| + |y|.
          VI is not referenced if WI = 0.0. \n
 * @param[out] B
          B is REAL array, dimension (LDB,N) \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= N+1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[in] EPS3
          EPS3 is REAL \n
          A small machine-dependent value which is used to perturb
          close eigenvalues, and to replace zero pivots. \n
 * @param[in] SMLNUM
          SMLNUM is REAL \n
          A machine-dependent value close to the underflow threshold. \n
 * @param[in] BIGNUM
          BIGNUM is REAL \n
          A machine-dependent value close to the overflow threshold. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          = 1:  inverse iteration did not converge; VR is set to the
                last iterate, and so is VI if WI.ne.0.0. \n

 * */
    template <typename T>
    void laein(logical *rightv, logical *noinit, integer *n, T *h, integer *ldh, T *wr, T *wi,
               T *vr, T *vi, T *b, integer *ldb, T *work, T *eps3, T *smlnum, T *bignum,
               integer *info)
    {
        laein(rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum, bignum, info);
    }
    template <typename T, typename Ta>
    void laein(logical *rightv, logical *noinit, integer *n, T *h, integer *ldh, T *w, T *v, T *b,
               integer *ldb, Ta *rwork, Ta *eps3, Ta *smlnum, integer *info)
    {
        laein(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info);
    }
    /** @}*/ // end of laein

    /** @defgroup lahqr lahqr
     * @ingroup NYS
     * @{
     */
    /*! @brief LAHQR computes the eigenvalues and Schur factorization of an upper \n
     Hessenberg matrix, using the double-shift/single-shift QR algorithm
 * @details
 * \b Purpose:
    \verbatim
    LAHQR is an auxiliary routine called by SHSEQR to update the
    eigenvalues and Schur decomposition already computed by SHSEQR, by
    dealing with the Hessenberg submatrix in rows and columns ILO to
    IHI.
    \endverbatim

 * @param[in] WANTT
          WANTT is LOGICAL \n
          = .TRUE. : the full Schur form T is required; \n
          = .FALSE.: only eigenvalues are required. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          = .TRUE. : the matrix of Schur vectors Z is required; \n
          = .FALSE.: Schur vectors are not required. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix H.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          It is assumed that H is already upper quasi-triangular in
          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
          ILO = 1). SLAHQR works primarily with the Hessenberg
          submatrix in rows and columns ILO to IHI, but applies
          transformations to all of H if WANTT is .TRUE.. \n
          1 <= ILO <= fla_max(1,IHI); IHI <= N. \n
 * @param[in,out] H
          H is REAL array, dimension (LDH,N) \n
          On entry, the upper Hessenberg matrix H. \n
          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
          quasi-triangular in rows and columns ILO:IHI, with any
          2-by-2 diagonal blocks in standard form. If INFO is zero
          and WANTT is .FALSE., the contents of H are unspecified on
          exit.  The output state of H if INFO is nonzero is given
          below under the description of INFO. \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the array H. LDH >= fla_max(1,N). \n
 * @param[out] WR
          WR is REAL array, dimension (N) \n
 * @param[out] WI
          WI is REAL array, dimension (N) \n
          The real and imaginary parts, respectively, of the computed
          eigenvalues ILO to IHI are stored in the corresponding
          elements of WR and WI. If two eigenvalues are computed as a
          complex conjugate pair, they are stored in consecutive
          elements of WR and WI, say the i-th and (i+1)th, with
          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
          eigenvalues are stored in the same order as on the diagonal
          of the Schur form   returned in H, with WR(i) = H(i,i), and, if
          H(i:i+1,i:i+1) is a 2-by-2 diagonal block, \n
          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i). \n
 * @param[in] ILOZ
          ILOZ is INTEGER \n
 * @param[in] IHIZ
          IHIZ is INTEGER \n
          Specify the rows of Z to which transformations must be
          applied if WANTZ is .TRUE.. \n
          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,N) \n
          If WANTZ is .TRUE., on entry Z must contain the current
          matrix Z of transformations accumulated by SHSEQR, and on
          exit Z has been updated; transformations are applied only to
          the submatrix Z(ILOZ:IHIZ,ILO:IHI). \n
          If WANTZ is .FALSE., Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z. LDZ >= fla_max(1,N). \n
 * @param[out] INFO
          INFO is INTEGER \n
           = 0:   successful exit \n
           > 0:   If INFO = i, SLAHQR failed to compute all the
                  eigenvalues ILO to IHI in a total of 30 iterations
                  per eigenvalue; elements i+1:ihi of WR and WI
                  contain those eigenvalues which have been
                  successfully computed. \n
 \n
                  If INFO > 0 and WANTT is .FALSE., then on exit,
                  the remaining unconverged eigenvalues are the
                  eigenvalues of the upper Hessenberg matrix rows
                  and columns ILO through INFO of the final, output
                  value of H. \n
 \n
                  If INFO > 0 and WANTT is .TRUE., then on exit
          (*)       (initial value of H)*U  = U*(final value of H)
                  where U is an orthogonal matrix.    The final
                  value of H is upper Hessenberg and triangular in
                  rows and columns INFO+1 through IHI. \n
 \n
                  If INFO > 0 and WANTZ is .TRUE., then on exit
                      (final value of Z)  = (initial value of Z)*U
                  where U is the orthogonal matrix in (*)
                  (regardless of the value of WANTT.)  \n

 * */
    template <typename T>
    void lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *wr, T *wi, integer *iloz, integer *ihiz, T *z, integer *ldz,
               integer *info)
    {
        lahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);
    }
    template <typename T>
    void lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *w, integer *iloz, integer *ihiz, T *z, integer *ldz, integer *info)
    {
        lahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
    }
    /** @}*/ // end of lahqr

    /** @defgroup laqr0 laqr0
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices
 from the Schur decomposition

 * @details
 * \b Purpose:
    \verbatim
    LAQR0 computes the eigenvalues of a Hessenberg matrix H
    and, optionally, the matrices T and Z from the Schur decomposition
    H = Z T Z**T, where T is an upper quasi-triangular matrix (the
    Schur form), and Z is the orthogonal matrix of Schur vectors.

    Optionally Z may be postmultiplied into an input orthogonal
    matrix Q so that this routine can give the Schur factorization
    of a matrix A which has been reduced to the Hessenberg form H
    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
    \endverbatim

 * @param[in] WANTT
          WANTT is LOGICAL \n
          = .TRUE. : the full Schur form T is required; \n
          = .FALSE.: only eigenvalues are required. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          = .TRUE. : the matrix of Schur vectors Z is required; \n
          = .FALSE.: Schur vectors are not required. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix H.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          It is assumed that H is already upper triangular in rows
          and columns 1:ILO-1 and IHI+1:N and, if ILO > 1,
          H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
          previous call to SGEBAL, and then passed to SGEHRD when the
          matrix output by SGEBAL is reduced to Hessenberg form.
          Otherwise, ILO and IHI should be set to 1 and N,
          respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
          If N = 0, then ILO = 1 and IHI = 0. \n
 * @param[in,out] H
          H is REAL array, dimension (LDH,N) \n
          On entry, the upper Hessenberg matrix H. \n
          On exit, if INFO = 0 and WANTT is .TRUE., then H contains
          the upper quasi-triangular matrix T from the Schur
          decomposition (the Schur form); 2-by-2 diagonal blocks
          (corresponding to complex conjugate pairs of eigenvalues)
          are   returned in standard form, with H(i,i) = H(i+1,i+1)
          and H(i+1,i)*H(i,i+1) < 0. If INFO = 0 and WANTT is
          .FALSE., then the contents of H are unspecified on exit.
          (The output value of H when INFO > 0 is given under the
          description of INFO below.) \n
 \n
          This subroutine may explicitly set H(i,j) = 0 for i > j and
          j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N. \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the array H. LDH >= fla_max(1,N). \n
 * @param[out] WR
          WR is REAL array, dimension (IHI) \n
 * @param[out] WI
          WI is REAL array, dimension (IHI) \n
          The real and imaginary parts, respectively, of the computed
          eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)
          and WI(ILO:IHI). If two eigenvalues are computed as a
          complex conjugate pair, they are stored in consecutive
          elements of WR and WI, say the i-th and (i+1)th, with
          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., then
          the eigenvalues are stored in the same order as on the
          diagonal of the Schur form   returned in H, with
          WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
          block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
          WI(i+1) = -WI(i). \n
 * @param[in] ILOZ
          ILOZ is INTEGER \n
 * @param[in] IHIZ
          IHIZ is INTEGER \n
          Specify the rows of Z to which transformations must be
          applied if WANTZ is .TRUE.. \n
          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,IHI) \n
          If WANTZ is .FALSE., then Z is not referenced. \n
          If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
          replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
          orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
          (The output value of Z when INFO > 0 is given under
          the description of INFO below.) \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  if WANTZ is .TRUE.
          then LDZ >= MAX(1,IHIZ).  Otherwise, LDZ >= 1. \n
 * @param[out] WORK
          WORK is REAL array, dimension LWORK \n
          On exit, if LWORK = -1, WORK(1)   returns an estimate of
          the optimal value for LWORK. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N)
          is sufficient, but LWORK typically as large as 6*N may
          be required for optimal performance.  A workspace query
          to determine the optimal workspace size is recommended.
 \n
          If LWORK = -1, then SLAQR0 does a workspace query.
          In this case, SLAQR0 checks the input parameters and
          estimates the optimal workspace size for the given
          values of N, ILO and IHI.  The estimate is   returned
          in WORK(1).  No error message related to LWORK is
          issued by XERBLA.  Neither H nor Z are accessed. \n
 * @param[out] INFO
          INFO is INTEGER
             = 0:  successful exit \n
             > 0:  if INFO = i, SLAQR0 failed to compute all of
                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
                and WI contain those eigenvalues which have been
                successfully computed.  (Failures are rare.)
 \n
                If INFO > 0 and WANT is .FALSE., then on exit,
                the remaining unconverged eigenvalues are the eigen-
                values of the upper Hessenberg matrix rows and
                columns ILO through INFO of the final, output
                value of H.
 \n
                If INFO > 0 and WANTT is .TRUE., then on exit
 \n
           (*)  (initial value of H)*U  = U*(final value of H)
 \n
                where U is an orthogonal matrix.  The final
                value of H is upper Hessenberg and quasi-triangular
                in rows and columns INFO+1 through IHI.
 \n
                If INFO > 0 and WANTZ is .TRUE., then on exit
 \n
                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
 \n
                where U is the orthogonal matrix in (*) (regard-
                less of the value of WANTT.)
 \n
                If INFO > 0 and WANTZ is .FALSE., then Z is not
                accessed. \n

 * */
    template <typename T>
    void laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *wr, T *wi, integer *iloz, integer *ihiz, T *z, integer *ldz,
               T *work, integer *lwork, integer *info)
    {
        laqr0(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
    }
    template <typename T>
    void laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *w, integer *iloz, integer *ihiz, T *z, integer *ldz, T *work,
               integer *lwork, integer *info)
    {
        laqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
    }
    /** @}*/ // end of laqr0

    /** @defgroup laqr1 laqr1
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3
 matrix H and specified shifts

 * @details
 * \b Purpose:
    \verbatim
     Given a 2-by-2 or 3-by-3 matrix H, SLAQR1 sets v to a
     scalar multiple of the first column of the product

     (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)

     scaling to avoid overflows and most underflows. It
     is assumed that either

             1) sr1 = sr2 and si1 = -si2
         or
             2) si1 = si2 = 0.

     This is useful for starting double implicit shift bulges
     in the QR algorithm.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          Order of the matrix H. N must be either 2 or 3. \n
 * @param[in] H
          H is REAL array, dimension (LDH,N) \n
          The 2-by-2 or 3-by-3 matrix H in (*). \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of H as declared in
          the calling procedure.  LDH >= N \n
 * @param[in] SR1
          SR1 is REAL \n
 * @param[in] SI1
          SI1 is REAL \n
 * @param[in] SR2
          SR2 is REAL \n
 * @param[in] SI2
          SI2 is REAL \n
          The shifts in (*). \n
 * @param[out] V
          V is REAL array, dimension (N) \n
          A scalar multiple of the first column of the
          matrix K in (*).  \n

 * */
    template <typename T>
    void laqr1(integer *n, T *h, integer *ldh, T *sr1, T *si1, T *sr2, T *si2, T *v)
    {
        laqr1(n, h, ldh, sr1, si1, sr2, si2, v);
    }
    template <typename T>
    void laqr1(integer *n, T *h, integer *ldh, T *s1, T *s2, T *v)
    {
        laqr1(n, h, ldh, s1, s2, v);
    }
    /** @}*/ // end of laqr1

    /** @defgroup laqr2 laqr2
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQR2 performs the orthogonal similarity transformation of a Hessenberg  \n
     matrix  to detect and deflate fully converged eigenvalues from a trailing      \n
     principal submatrix (aggressive early deflation)
 * @details
 * \b Purpose:
    \verbatim
    LAQR2 is identical to SLAQR3 except that it avoids
    recursion by calling SLAHQR instead of SLAQR4.

    Aggressive early deflation:

    This subroutine accepts as input an upper Hessenberg matrix
    H and performs an orthogonal similarity transformation
    designed to detect and deflate fully converged eigenvalues from
    a trailing principal submatrix.  On output H has been over-
    written by a new Hessenberg matrix that is a perturbation of
    an orthogonal similarity transformation of H.  It is to be
    hoped that the final version of H has many zero subdiagonal
    entries.
    \endverbatim

 * @param[in] WANTT
          WANTT is LOGICAL \n
          If .TRUE., then the Hessenberg matrix H is fully updated
          so that the quasi-triangular Schur factor may be
          computed (in cooperation with the calling subroutine).
          If .FALSE., then only enough of H is updated to preserve
          the eigenvalues. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          If .TRUE., then the orthogonal matrix Z is updated so
          so that the orthogonal Schur factor may be computed
          (in cooperation with the calling subroutine).
          If .FALSE., then Z is not referenced. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix H and (if WANTZ is .TRUE.) the
          order of the orthogonal matrix Z. \n
 * @param[in] KTOP
          KTOP is INTEGER \n
          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
          KBOT and KTOP together determine an isolated block
          along the diagonal of the Hessenberg matrix. \n
 * @param[in] KBOT
          KBOT is INTEGER \n
          It is assumed without a check that either
          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
          determine an isolated block along the diagonal of the
          Hessenberg matrix. \n
 * @param[in] NW
          NW is INTEGER \n
          Deflation window size.  1 <= NW <= (KBOT-KTOP+1). \n
 * @param[in,out] H
          H is REAL array, dimension (LDH,N) \n
          On input the initial N-by-N section of H stores the
          Hessenberg matrix undergoing aggressive early deflation. \n
          On output H has been transformed by an orthogonal
          similarity transformation, perturbed, and the   returned
          to Hessenberg form that (it is to be hoped) has some
          zero subdiagonal entries. \n
 * @param[in] LDH
          LDH is INTEGER \n
          Leading dimension of H just as declared in the calling
          subroutine.  N <= LDH \n
 * @param[in] ILOZ
          ILOZ is INTEGER \n
 * @param[in] IHIZ
          IHIZ is INTEGER \n
          Specify the rows of Z to which transformations must be
          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,N) \n
          IF WANTZ is .TRUE., then on output, the orthogonal
          similarity transformation mentioned above has been
          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
          If WANTZ is .FALSE., then Z is unreferenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of Z just as declared in the
          calling subroutine.  1 <= LDZ. \n
 * @param[out] NS
          NS is INTEGER \n
          The number of unconverged (ie approximate) eigenvalues
          returned in SR and SI that may be used as shifts by the
          calling subroutine. \n
 * @param[out] ND
          ND is INTEGER \n
          The number of converged eigenvalues uncovered by this
          subroutine. \n
 * @param[out] SR
          SR is REAL array, dimension (KBOT) \n
 * @param[out] SI
          SI is REAL array, dimension (KBOT) \n
          On output, the real and imaginary parts of approximate
          eigenvalues that may be used for shifts are stored in
          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
          The real and imaginary parts of converged eigenvalues
          are stored in SR(KBOT-ND+1) through SR(KBOT) and
          SI(KBOT-ND+1) through SI(KBOT), respectively. \n
 * @param[out] V
          V is REAL array, dimension (LDV,NW) \n
          An NW-by-NW work array. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of V just as declared in the
          calling subroutine.  NW <= LDV \n
 * @param[in] NH
          NH is INTEGER \n
          The number of columns of T.  NH >= NW. \n
 * @param[out] T
          T is REAL array, dimension (LDT,NW) \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of T just as declared in the
          calling subroutine.  NW <= LDT \n
 * @param[in] NV
          NV is INTEGER \n
          The number of rows of work array WV available for
          workspace.  NV >= NW. \n
 * @param[out] WV
          WV is REAL array, dimension (LDWV,NW) \n
 * @param[in] LDWV
          LDWV is INTEGER \n
          The leading dimension of W just as declared in the
          calling subroutine.  NW <= LDV \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
          On exit, WORK(1) is set to an estimate of the optimal value
          of LWORK for the given values of N, NW, KTOP and KBOT. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the work array WORK.  LWORK = 2*NW
          suffices, but greater efficiency may result from larger
          values of LWORK. \n
 \n
          If LWORK = -1, then a workspace query is assumed; SLAQR2
          only estimates the optimal workspace size for the given
          values of N, NW, KTOP and KBOT.  The estimate is   returned
          in WORK(1).  No error message related to LWORK is issued
          by XERBLA.  Neither H nor Z are accessed. \n

 * */
    template <typename T>
    void laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
               integer *nw, T *h, integer *ldh, integer *iloz, integer *ihiz, T *z, integer *ldz,
               integer *ns, integer *nd, T *sr, T *si, T *v, integer *ldv, integer *nh, T *t,
               integer *ldt, integer *nv, T *wv, integer *ldwv, T *work, integer *lwork)
    {
        laqr2(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv,
              nh, t, ldt, nv, wv, ldwv, work, lwork);
    }
    template <typename T>
    void laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
               integer *nw, T *h, integer *ldh, integer *iloz, integer *ihiz, T *z, integer *ldz,
               integer *ns, integer *nd, T *sh, T *v, integer *ldv, integer *nh, T *t, integer *ldt,
               integer *nv, T *wv, integer *ldwv, T *work, integer *lwork)
    {
        laqr2(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh,
              t, ldt, nv, wv, ldwv, work, lwork);
    }
    /** @}*/ // end of laqr2

    /** @defgroup laqr3 laqr3
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQR3 performs the orthogonal similarity transformation of a Hessenberg \n
     matrix to detect and deflate fully converged eigenvalues from a trailing      \n
     principal submatrix (aggressive early deflation).
 * @details
 * \b Purpose:
    \verbatim
    LAQR3 accepts as input an upper Hessenberg matrix
    H and performs an orthogonal similarity transformation
    designed to detect and deflate fully converged eigenvalues from
    a trailing principal submatrix.  On output H has been over-
    written by a new Hessenberg matrix that is a perturbation of
    an orthogonal similarity transformation of H.  It is to be
    hoped that the final version of H has many zero subdiagonal
    entries.
    \endverbatim

 * @param[in] WANTT
          WANTT is LOGICAL \n
          If .TRUE., then the Hessenberg matrix H is fully updated
          so that the quasi-triangular Schur factor may be
          computed (in cooperation with the calling subroutine). \n
          If .FALSE., then only enough of H is updated to preserve
          the eigenvalues. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          If .TRUE., then the orthogonal matrix Z is updated so
          so that the orthogonal Schur factor may be computed
          (in cooperation with the calling subroutine). \n
          If .FALSE., then Z is not referenced. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix H and (if WANTZ is .TRUE.) the
          order of the orthogonal matrix Z. \n
 * @param[in] KTOP
          KTOP is INTEGER \n
          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
          KBOT and KTOP together determine an isolated block
          along the diagonal of the Hessenberg matrix. \n
 * @param[in] KBOT
          KBOT is INTEGER \n
          It is assumed without a check that either
          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
          determine an isolated block along the diagonal of the
          Hessenberg matrix. \n
 * @param[in] NW
          NW is INTEGER \n
          Deflation window size.  1 <= NW <= (KBOT-KTOP+1). \n
 * @param[in,out] H
          H is REAL array, dimension (LDH,N) \n
          On input the initial N-by-N section of H stores the
          Hessenberg matrix undergoing aggressive early deflation.
          On output H has been transformed by an orthogonal
          similarity transformation, perturbed, and the   returned
          to Hessenberg form that (it is to be hoped) has some
          zero subdiagonal entries. \n
 * @param[in] LDH
          LDH is INTEGER \n
          Leading dimension of H just as declared in the calling
          subroutine.  N <= LDH \n
 * @param[in] ILOZ
          ILOZ is INTEGER \n
 * @param[in] IHIZ
          IHIZ is INTEGER \n
          Specify the rows of Z to which transformations must be
          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,N) \n
          IF WANTZ is .TRUE., then on output, the orthogonal
          similarity transformation mentioned above has been
          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right. \n
          If WANTZ is .FALSE., then Z is unreferenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of Z just as declared in the
          calling subroutine.  1 <= LDZ. \n
 * @param[out] NS
          NS is INTEGER \n
          The number of unconverged (ie approximate) eigenvalues
          returned in SR and SI that may be used as shifts by the
          calling subroutine. \n
 * @param[out] ND
          ND is INTEGER \n
          The number of converged eigenvalues uncovered by this
          subroutine. \n
 * @param[out] SR
          SR is REAL array, dimension (KBOT) \n
 * @param[out] SI
          SI is REAL array, dimension (KBOT) \n
          On output, the real and imaginary parts of approximate
          eigenvalues that may be used for shifts are stored in
          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
          The real and imaginary parts of converged eigenvalues
          are stored in SR(KBOT-ND+1) through SR(KBOT) and
          SI(KBOT-ND+1) through SI(KBOT), respectively. \n
 * @param[out] V
          V is REAL array, dimension (LDV,NW) \n
          An NW-by-NW work array. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of V just as declared in the
          calling subroutine.  NW <= LDV \n
 * @param[in] NH
          NH is INTEGER \n
          The number of columns of T.  NH >= NW. \n
 * @param[out] T
          T is REAL array, dimension (LDT,NW) \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of T just as declared in the
          calling subroutine.  NW <= LDT \n
 * @param[in] NV
          NV is INTEGER \n
          The number of rows of work array WV available for
          workspace.  NV >= NW. \n
 * @param[out] WV
          WV is REAL array, dimension (LDWV,NW) \n
 * @param[in] LDWV
          LDWV is INTEGER \n
          The leading dimension of W just as declared in the
          calling subroutine.  NW <= LDV \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
          On exit, WORK(1) is set to an estimate of the optimal value
          of LWORK for the given values of N, NW, KTOP and KBOT. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the work array WORK.  LWORK = 2*NW
          suffices, but greater efficiency may result from larger
          values of LWORK. \n
 \n
          If LWORK = -1, then a workspace query is assumed; SLAQR3
          only estimates the optimal workspace size for the given
          values of N, NW, KTOP and KBOT.  The estimate is   returned
          in WORK(1).  No error message related to LWORK is issued
          by XERBLA.  Neither H nor Z are accessed.  \n

 * */
    template <typename T>
    void laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
               integer *nw, T *h, integer *ldh, integer *iloz, integer *ihiz, T *z, integer *ldz,
               integer *ns, integer *nd, T *sr, T *si, T *v, integer *ldv, integer *nh, T *t,
               integer *ldt, integer *nv, T *wv, integer *ldwv, T *work, integer *lwork)
    {
        laqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv,
              nh, t, ldt, nv, wv, ldwv, work, lwork);
    }
    template <typename T>
    void laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
               integer *nw, T *h, integer *ldh, integer *iloz, integer *ihiz, T *z, integer *ldz,
               integer *ns, integer *nd, T *sh, T *v, integer *ldv, integer *nh, T *t, integer *ldt,
               integer *nv, T *wv, integer *ldwv, T *work, integer *lwork)
    {
        laqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh,
              t, ldt, nv, wv, ldwv, work, lwork);
    }
    /** @}*/ // end of laqr3

    /** @defgroup laqr4 laqr4
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices
 from the Schur decomposition

 * @details
 * \b Purpose:
    \verbatim
    LAQR4 implements one level of recursion for SLAQR0.
    It is a complete implementation of the small bulge multi-shift
    QR algorithm.  It may be called by SLAQR0 and, for large enough
    deflation window size, it may be called by SLAQR3.  This
    subroutine is identical to SLAQR0 except that it calls SLAQR2
    instead of SLAQR3.

    SLAQR4 computes the eigenvalues of a Hessenberg matrix H
    and, optionally, the matrices T and Z from the Schur decomposition
    H = Z T Z**T, where T is an upper quasi-triangular matrix (the
    Schur form), and Z is the orthogonal matrix of Schur vectors.

    Optionally Z may be postmultiplied into an input orthogonal
    matrix Q so that this routine can give the Schur factorization
    of a matrix A which has been reduced to the Hessenberg form H
    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
    \endverbatim

 * @param[in] WANTT
          WANTT is LOGICAL \n
          = .TRUE. : the full Schur form T is required; \n
          = .FALSE.: only eigenvalues are required. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          = .TRUE. : the matrix of Schur vectors Z is required; \n
          = .FALSE.: Schur vectors are not required. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix H.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          It is assumed that H is already upper triangular in rows
          and columns 1:ILO-1 and IHI+1:N and, if ILO > 1,
          H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
          previous call to SGEBAL, and then passed to SGEHRD when the
          matrix output by SGEBAL is reduced to Hessenberg form.
          Otherwise, ILO and IHI should be set to 1 and N,
          respectively.  If N > 0, then 1 <= ILO <= IHI <= N. \n
          If N = 0, then ILO = 1 and IHI = 0. \n
 * @param[in,out] H
          H is REAL array, dimension (LDH,N) \n
          On entry, the upper Hessenberg matrix H. \n
          On exit, if INFO = 0 and WANTT is .TRUE., then H contains
          the upper quasi-triangular matrix T from the Schur
          decomposition (the Schur form); 2-by-2 diagonal blocks
          (corresponding to complex conjugate pairs of eigenvalues)
          are   returned in standard form, with H(i,i) = H(i+1,i+1)
          and H(i+1,i)*H(i,i+1) < 0. If INFO = 0 and WANTT is
          .FALSE., then the contents of H are unspecified on exit.
          (The output value of H when INFO > 0 is given under the
          description of INFO below.) \n
 \n
          This subroutine may explicitly set H(i,j) = 0 for i > j and
          j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N. \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the array H. LDH >= fla_max(1,N). \n
 * @param[out] WR
          WR is REAL array, dimension (IHI) \n
 * @param[out] WI
          WI is REAL array, dimension (IHI) \n
          The real and imaginary parts, respectively, of the computed
          eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)
          and WI(ILO:IHI). If two eigenvalues are computed as a
          complex conjugate pair, they are stored in consecutive
          elements of WR and WI, say the i-th and (i+1)th, with
          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., then
          the eigenvalues are stored in the same order as on the
          diagonal of the Schur form   returned in H, with
          WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
          block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
          WI(i+1) = -WI(i). \n
 * @param[in] ILOZ
          ILOZ is INTEGER \n
 * @param[in] IHIZ
          IHIZ is INTEGER \n
          Specify the rows of Z to which transformations must be
          applied if WANTZ is .TRUE.. \n
          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,IHI) \n
          If WANTZ is .FALSE., then Z is not referenced. \n
          If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
          replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
          orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
          (The output value of Z when INFO > 0 is given under
          the description of INFO below.) \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  if WANTZ is .TRUE.
          then LDZ >= MAX(1,IHIZ).  Otherwise, LDZ >= 1. \n
 * @param[out] WORK
          WORK is REAL array, dimension LWORK \n
          On exit, if LWORK = -1, WORK(1)   returns an estimate of
          the optimal value for LWORK. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N)
          is sufficient, but LWORK typically as large as 6*N may
          be required for optimal performance.  A workspace query
          to determine the optimal workspace size is recommended.
 \n
          If LWORK = -1, then SLAQR4 does a workspace query.
          In this case, SLAQR4 checks the input parameters and
          estimates the optimal workspace size for the given
          values of N, ILO and IHI.  The estimate is   returned
          in WORK(1).  No error message related to LWORK is
          issued by XERBLA.  Neither H nor Z are accessed. \n
 * @param[out] INFO
          INFO is INTEGER \n
             = 0:  successful exit \n
             > 0:  if INFO = i, SLAQR4 failed to compute all of
                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
                and WI contain those eigenvalues which have been
                successfully computed.  (Failures are rare.)
 \n
                If INFO > 0 and WANT is .FALSE., then on exit,
                the remaining unconverged eigenvalues are the eigen-
                values of the upper Hessenberg matrix rows and
                columns ILO through INFO of the final, output
                value of H.
 \n
                If INFO > 0 and WANTT is .TRUE., then on exit
 \n
           (*)  (initial value of H)*U  = U*(final value of H)
 \n
                where U is a orthogonal matrix.  The final
                value of  H is upper Hessenberg and triangular in
                rows and columns INFO+1 through IHI.
 \n
                If INFO > 0 and WANTZ is .TRUE., then on exit
 \n
                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
 \n
                where U is the orthogonal matrix in (*) (regard-
                less of the value of WANTT.)
 \n
                If INFO > 0 and WANTZ is .FALSE., then Z is not
                accessed.  \n

 * */
    template <typename T>
    void laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *wr, T *wi, integer *iloz, integer *ihiz, T *z, integer *ldz,
               T *work, integer *lwork, integer *info)
    {
        laqr4(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
    }
    template <typename T>
    void laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *w, integer *iloz, integer *ihiz, T *z, integer *ldz, T *work,
               integer *lwork, integer *info)
    {
        laqr4(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
    }
    /** @}*/ // end of laqr4

    /** @defgroup laqr5 laqr5
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQR5 performs a single small-bulge multi-shift QR sweep

 * @details
 * \b Purpose:
    \verbatim
    LAQR5, called by SLAQR0, performs a
    single small-bulge multi-shift QR sweep.
    \endverbatim

 * @param[in] WANTT
          WANTT is LOGICAL \n
          WANTT = .true. if the quasi-triangular Schur factor
          is being computed.  WANTT is set to .false. otherwise. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          WANTZ = .true. if the orthogonal Schur factor is being
          computed.  WANTZ is set to .false. otherwise. \n
 * @param[in] KACC22
          KACC22 is INTEGER with value 0, 1, or 2. \n
          Specifies the computation mode of far-from-diagonal
          orthogonal updates. \n
          = 0: SLAQR5 does not accumulate reflections and does not
               use matrix-matrix multiply to update far-from-diagonal
               matrix entries. \n
          = 1: SLAQR5 accumulates reflections and uses matrix-matrix
               multiply to update the far-from-diagonal matrix entries. \n
          = 2: SLAQR5 accumulates reflections, uses matrix-matrix
               multiply to update the far-from-diagonal matrix entries,
               and takes advantage of 2-by-2 block structure during
               matrix multiplies. \n
 * @param[in] N
          N is INTEGER \n
          N is the order of the Hessenberg matrix H upon which this
          subroutine operates. \n
 * @param[in] KTOP
          KTOP is INTEGER \n
 * @param[in] KBOT
          KBOT is INTEGER \n
          These are the first and last rows and columns of an
          isolated diagonal block upon which the QR sweep is to be
          applied. It is assumed without a check that \n
                   either KTOP = 1  or   H(KTOP,KTOP-1) = 0 \n
          and \n
                   either KBOT = N  or   H(KBOT+1,KBOT) = 0. \n
 * @param[in] NSHFTS
          NSHFTS is INTEGER \n
          NSHFTS gives the number of simultaneous shifts.  NSHFTS
          must be positive and even. \n
 * @param[in,out] SR
          SR is REAL array, dimension (NSHFTS) \n
 * @param[in,out] SI
          SI is REAL array, dimension (NSHFTS) \n
          SR contains the real parts and SI contains the imaginary
          parts of the NSHFTS shifts of origin that define the
          multi-shift QR sweep.  On output SR and SI may be
          reordered. \n
 * @param[in,out] H
          H is REAL array, dimension (LDH,N) \n
          On input H contains a Hessenberg matrix.  On output a
          multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
          to the isolated diagonal block in rows and columns KTOP
          through KBOT. \n
 * @param[in] LDH
          LDH is INTEGER \n
          LDH is the leading dimension of H just as declared in the
          calling procedure.  LDH >= MAX(1,N). \n
 * @param[in] ILOZ
          ILOZ is INTEGER \n
 * @param[in] IHIZ
          IHIZ is INTEGER \n
          Specify the rows of Z to which transformations must be
          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,IHIZ) \n
          If WANTZ = .TRUE., then the QR Sweep orthogonal
          similarity transformation is accumulated into
          Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right. \n
          If WANTZ = .FALSE., then Z is unreferenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          LDZ is the leading dimension of Z just as declared in
          the calling procedure. LDZ >= N. \n
 * @param[out] V
          V is REAL array, dimension (LDV,NSHFTS/2) \n
 * @param[in] LDV
          LDV is INTEGER \n
          LDV is the leading dimension of V as declared in the
          calling procedure.  LDV >= 3. \n
 * @param[out] U
          U is REAL array, dimension (LDU,3*NSHFTS-3) \n
 * @param[in] LDU
          LDU is INTEGER \n
          LDU is the leading dimension of U just as declared in the
          in the calling subroutine.  LDU >= 3*NSHFTS-3. \n
 * @param[in] NV
          NV is INTEGER \n
          NV is the number of rows in WV agailable for workspace.
          NV >= 1. \n
 * @param[out] WV
          WV is REAL array, dimension (LDWV,3*NSHFTS-3)
 * @param[in] LDWV
          LDWV is INTEGER
          LDWV is the leading dimension of WV as declared in the
          in the calling subroutine.  LDWV >= NV.
 * @param[in] NH
          NH is INTEGER
          NH is the number of columns in array WH available for
          workspace. NH >= 1.
 * @param[out] WH
          WH is REAL array, dimension (LDWH,NH) \n
 * @param[in] LDWH
          LDWH is INTEGER \n
          Leading dimension of WH just as declared in the
          calling procedure.  LDWH >= 3*NSHFTS-3.  \n

 * */
    template <typename T>
    void laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
               integer *kbot, integer *nshfts, T *sr, T *si, T *h, integer *ldh, integer *iloz,
               integer *ihiz, T *z, integer *ldz, T *v, integer *ldv, T *u, integer *ldu,
               integer *nv, T *wv, integer *ldwv, integer *nh, T *wh, integer *ldwh)
    {
        laqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v,
              ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
    }
    template <typename T>
    void laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
               integer *kbot, integer *nshfts, T *s, T *h, integer *ldh, integer *iloz,
               integer *ihiz, T *z, integer *ldz, T *v, integer *ldv, T *u, integer *ldu,
               integer *nv, T *wv, integer *ldwv, integer *nh, T *wh, integer *ldwh)
    {
        laqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u,
              ldu, nv, wv, ldwv, nh, wh, ldwh);
    }
    /** @}*/ // end of laqr5

    /** @defgroup laqz0 laqz0
     * @ingroup NYS
     * @{
     */
    /*! @brief  Computes the eigenvalues of a real matrix pair (H,T),
            where H is an upper Hessenberg matrix and T is upper triangular,
            using the double-shift QZ method.

 * @details
 * \b Purpose:
    \verbatim
    LAQZ0 computes the eigenvalues of a real matrix pair (H,T),
    where H is an upper Hessenberg matrix and T is upper triangular,
    using the double-shift QZ method.
    Matrix pairs of this type are produced by the reduction to
    generalized upper Hessenberg form of a real matrix pair (A,B):

      A = Q1*H*Z1**T,  B = Q1*T*Z1**T,

    as computed by SGGHRD.

    If JOB='S', then the Hessenberg-triangular pair (H,T) is
    also reduced to generalized Schur form,

      H = Q*S*Z**T,  T = Q*P*Z**T,

    where Q and Z are orthogonal matrices, P is an upper triangular
    matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2
    diagonal blocks.

    The 1-by-1 blocks correspond to real eigenvalues of the matrix pair
    (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of
    eigenvalues.

    Additionally, the 2-by-2 upper triangular diagonal blocks of P
    corresponding to 2-by-2 blocks of S are reduced to positive diagonal
    form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,
    P(j,j) > 0, and P(j+1,j+1) > 0.

    Optionally, the orthogonal matrix Q from the generalized Schur
    factorization may be postmultiplied into an input matrix Q1, and the
    orthogonal matrix Z may be postmultiplied into an input matrix Z1.
    If Q1 and Z1 are the orthogonal matrices from SGGHRD that reduced
    the matrix pair (A,B) to generalized upper Hessenberg form, then the
    output matrices Q1*Q and Z1*Z are the orthogonal factors from the
    generalized Schur factorization of (A,B):

      A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.

    To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,
    of (A,B)) are computed as a pair of values (alpha,beta), where alpha is
    complex and beta real.
    If beta is nonzero, lambda = alpha / beta is an eigenvalue of the
    generalized nonsymmetric eigenvalue problem (GNEP)
      A*x = lambda*B*x
    and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
    alternate form of the GNEP
      mu*A*y = B*y.
    Real eigenvalues can be read directly from the generalized Schur
    form:
     alpha = S(i,i), beta = P(i,i).

    Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
        Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
        pp. 241--256.

    Ref: B. Kagstrom, D. Kressner, "Multishift Variants of the QZ
        Algorithm with Aggressive Early Deflation", SIAM J. Numer.
        Anal., 29(2006), pp. 199--227.

    Ref: T. Steel, D. Camps, K. Meerbergen, R. Vandebril "A multishift,
        multipole rational QZ method with agressive early deflation"
    \endverbatim

 * @param[in] WANTS
          WANTS is CHARACTER*1 \n
          = 'E': Compute eigenvalues only; \n
          = 'S': Compute eigenvalues and the Schur form. \n
 * @param[in] WANTQ
          WANTQ is CHARACTER*1 \n
          = 'N': Left Schur vectors (Q) are not computed; \n
          = 'I': Q is initialized to the unit matrix and the matrix Q
                 of left Schur vectors of (A,B) is returned; \n
          = 'V': Q must contain an orthogonal matrix Q1 on entry and
                 the product Q1*Q is returned. \n
 * @param[in] WANTZ
          WANTZ is CHARACTER*1 \n
          = 'N': Right Schur vectors (Z) are not computed; \n
          = 'I': Z is initialized to the unit matrix and the matrix Z
                 of right Schur vectors of (A,B) is returned; \n
          = 'V': Z must contain an orthogonal matrix Z1 on entry and
                 the product Z1*Z is returned. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, Q, and Z.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          ILO and IHI mark the rows and columns of A which are in
          Hessenberg form.  It is assumed that A is already upper
          triangular in rows and columns 1:ILO-1 and IHI+1:N. \n
          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0. \n
 * @param[in, out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the N-by-N upper Hessenberg matrix A. \n
          On exit, if JOB = 'S', A contains the upper quasi-triangular
          matrix S from the generalized Schur factorization. \n
          If JOB = 'E', the diagonal blocks of A match those of S, but
          the rest of A is unspecified. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max( 1, N ). \n
 * @param[in, out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the N-by-N upper triangular matrix B. \n
          On exit, if JOB = 'S', B contains the upper triangular
          matrix P from the generalized Schur factorization; \n
          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S
          are reduced to positive diagonal form, i.e., if A(j+1,j) is
          non-zero, then B(j+1,j) = B(j,j+1) = 0, B(j,j) > 0, and
          B(j+1,j+1) > 0. \n
          If JOB = 'E', the diagonal blocks of B match those of P, but
          the rest of B is unspecified. \n
 * @param[in]	LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max( 1, N ). \n
 * @param[out]	ALPHAR
          ALPHAR is REAL array, dimension (N) \n
          The real parts of each scalar alpha defining an eigenvalue
          of GNEP. \n
 * @param[out]	ALPHAI
          ALPHAI is REAL array, dimension (N) \n
          The imaginary parts of each scalar alpha defining an
          eigenvalue of GNEP. \n
          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
          positive, then the j-th and (j+1)-st eigenvalues are a
          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j). \n
 * @param[out]	BETA
          BETA is REAL array, dimension (N) \n
          The scalars beta that define the eigenvalues of GNEP. \n
          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
          beta = BETA(j) represent the j-th eigenvalue of the matrix
          pair (A,B), in one of the forms lambda = alpha/beta or
          mu = beta/alpha.  Since either lambda or mu may overflow,
          they should not, in general, be computed. \n
 * @param[in, out] Q
          Q is REAL array, dimension (LDQ, N) \n
          On entry, if COMPQ = 'V', the orthogonal matrix Q1 used in
          the reduction of (A,B) to generalized Hessenberg form. \n
          On exit, if COMPQ = 'I', the orthogonal matrix of left Schur
          vectors of (A,B), and if COMPQ = 'V', the orthogonal matrix
          of left Schur vectors of (A,B). \n
          Not referenced if COMPQ = 'N'. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= 1. \n
          If COMPQ='V' or 'I', then LDQ >= N. \n
 * @param[in, out] Z
          Z is REAL array, dimension (LDZ, N) \n
          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in
          the reduction of (A,B) to generalized Hessenberg form. \n
          On exit, if COMPZ = 'I', the orthogonal matrix of
          right Schur vectors of (H,T), and if COMPZ = 'V', the
          orthogonal matrix of right Schur vectors of (A,B). \n
          Not referenced if COMPZ = 'N'. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1. \n
          If COMPZ='V' or 'I', then LDZ >= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N). \n

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[in] REC
          REC is INTEGER \n
          REC indicates the current recursion level. Should be set
          to 0 on first call. \n
 * @param[in] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void laqz0(char *wants, char *wantq, char *wantz, integer *n, integer *ilo, integer *ihi, T *a,
               integer *lda, T *b, integer *ldb, T *alphar, T *alphai, T *beta, T *q, integer *ldq,
               T *z, integer *ldz, T *work, integer *lwork, integer *rec, integer *info)
    {
        laqz0(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z,
              ldz, work, lwork, rec, info);
    }
    template <typename T, typename Ta>
    void laqz0(char *wants, char *wantq, char *wantz, integer *n, integer *ilo, integer *ihi, T *a,
               integer *lda, T *b, integer *ldb, T *alpha, T *beta, T *q, integer *ldq, T *z,
               integer *ldz, T *work, integer *lwork, Ta *rwork, integer *rec, integer *info)
    {
        laqz0(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, work,
              lwork, rwork, rec, info);
    }
    /** @}*/ // end of laqz0

    /** @defgroup laqz1 laqz1
     * @ingroup NYS
     * @{
     */
    /*! @brief Given a 3-by-3 matrix pencil (A,B), LAQZ1 sets v to a
           scalar multiple of the first column of the product.

 * @details
 * \b Purpose:
    \verbatim
    Given a 3-by-3 matrix pencil (A,B), LAQZ1 sets v to a
      scalar multiple of the first column of the product

      (*)  K = (A - (beta2*sr2 - i*si)*B)*B^(-1)*(beta1*A - (sr2 + i*si2)*B)*B^(-1).

      It is assumed that either

              1) sr1 = sr2
          or
              2) si = 0.

      This is useful for starting double implicit shift bulges
      in the QZ algorithm.
    \endverbatim

 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The 3-by-3 matrix A in (*). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A as declared in
          the calling procedure. \n
 * @param[in] B
          B is REAL array, dimension (LDB,N) \n
          The 3-by-3 matrix B in (*). \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B as declared in
          the calling procedure. \n
 * @param[in] SR1
          SR1 is REAL \n
 * @param[in] SR2
          SR2 is REAL \n
 * @param[in] SI
          SI is REAL \n
 * @param[in] BETA1
          BETA1 is REAL \n
 * @param[in] BETA2
          BETA2 is REAL \n
 * @param[out] V
          V is REAL array, dimension (N) \n
          A scalar multiple of the first column of the
          matrix K in (*). \n

 * */
    template <typename T>
    void laqz1(T *a, integer *lda, T *b, integer *ldb, T *sr1, T *sr2, T *si, T *beta1, T *beta2,
               T *v)
    {
        laqz1(a, lda, b, ldb, sr1, sr2, si, beta1, beta2, v);
    }
    template <typename T>
    void laqz1(logical *ilq, logical *ilz, integer *k, integer *istartm, integer *istopm,
               integer *ihi, T *a, integer *lda, T *b, integer *ldb, integer *nq, integer *qstart,
               T *q, integer *ldq, integer *nz, integer *zstart, T *z, integer *ldz)
    {
        laqz1(ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz, zstart, z,
              ldz);
    }
    /** @}*/ // end of laqz1

    /** @defgroup laqz2 laqz2
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQZ2 chases a 2x2 shift bulge in a matrix pencil down a single position.

 * @details
 * \b Purpose:
    \verbatim
    LAQZ2 chases a 2x2 shift bulge in a matrix pencil down a single position
    \endverbatim

 * @param[in] ILQ
          ILQ is LOGICAL \n
          Determines whether or not to update the matrix Q \n
 * @param[in] ILZ
          ILZ is LOGICAL \n
          Determines whether or not to update the matrix Z \n
 * @param[in] K
          K is INTEGER \n
          Index indicating the position of the bulge. \n
          On entry, the bulge is located in
          (A(k+1:k+2,k:k+1),B(k+1:k+2,k:k+1)). \n
          On exit, the bulge is located in
          (A(k+2:k+3,k+1:k+2),B(k+2:k+3,k+1:k+2)). \n
 * @param[in] ISTARTM
          ISTARTM is INTEGER \n
 * @param[in] ISTOPM
          ISTOPM is INTEGER \n
          Updates to (A,B) are restricted to
          (istartm:k+3,k:istopm). It is assumed
          without checking that istartm <= k+1 and
          k+2 <= istopm \n
 * @param[in] IHI
          IHI is INTEGER \n
 * @param[in, out] A
          A is REAL array, dimension (LDA,N) \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of A as declared in
          the calling procedure. \n
 * @param[in, out] B
          B is REAL array, dimension (LDB,N) \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B as declared in
          the calling procedure. \n
 * @param[in] NQ
          NQ is INTEGER \n
          The order of the matrix Q \n
 * @param[in] QSTART
          QSTART is INTEGER \n
          Start index of the matrix Q. Rotations are applied
          To columns k+2-qStart:k+4-qStart of Q. \n
 * @param[in, out] Q
          Q is REAL array, dimension (LDQ,NQ) \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of Q as declared in
          the calling procedure. \n
 * @param[in] NZ
          NZ is INTEGER \n
          The order of the matrix Z \n
 * @param[in] ZSTART
          ZSTART is INTEGER \n
          Start index of the matrix Z. Rotations are applied
          To columns k+1-qStart:k+3-qStart of Z. \n
 * @param[in, out] Z
          Z is REAL array, dimension (LDZ,NZ) \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of Q as declared in
          the calling procedure. \n

 * */
    template <typename T>
    void laqz2(logical *ilq, logical *ilz, integer *k, integer *istartm, integer *istopm,
               integer *ihi, T *a, integer *lda, T *b, integer *ldb, integer *nq, integer *qstart,
               T *q, integer *ldq, integer *nz, integer *zstart, T *z, integer *ldz)
    {
        laqz2(ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz, zstart, z,
              ldz);
    }
    template <typename T, typename Ta>
    void laqz2(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi,
               integer *nw, T *a, integer *lda, T *b, integer *ldb, T *q, integer *ldq, T *z,
               integer *ldz, integer *ns, integer *nd, T *alpha, T *beta, T *qc, integer *ldqc,
               T *zc, integer *ldzc, T *work, integer *lwork, float *rwork, integer *rec,
               integer *info)
    {
        laqz2(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, ns, nd, alpha,
              beta, qc, ldqc, zc, ldzc, work, lwork, rwork, rec, info);
    }
    /** @}*/ // end of laqz2

    /** @defgroup laqz3 laqz3
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQZ3 performs AED

 * @details
 * \b Purpose:
    \verbatim
    LAQZ3 performs AED
    \endverbatim

 * @param[in] ILSCHUR
          ILSCHUR is LOGICAL \n
          Determines whether or not to update the full Schur form \n
 * @param[in] ILQ
          ILQ is LOGICAL \n
          Determines whether or not to update the matrix Q \n
 * @param[in] ILZ
          ILZ is LOGICAL \n
          Determines whether or not to update the matrix Z \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, Q, and Z.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          ILO and IHI mark the rows and columns of (A,B) which
          are to be normalized \n
 * @param[in] NW
          NW is INTEGER \n
          The desired size of the deflation window. \n
 * @param[in, out] A
          A is REAL array, dimension (LDA, N) \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max( 1, N ). \n
 * @param[in, out] B
          B is REAL array, dimension (LDB,N) \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max( 1, N ). \n
 * @param[in, out] Q
          Q is REAL array, dimension (LDQ,N) \n
 * @param[in] LDQ
          LDQ is INTEGER \n
 * @param[in, out] Z
          Z is REAL array, dimension (LDZ,N) \n
 * @param[in] LDZ
          LDZ is INTEGER \n
 * @param[out] NS
          NS is INTEGER \n
          The number of unconverged eigenvalues available to
          use as shifts. \n
 * @param[out] ND
          ND is INTEGER \n
          The number of converged eigenvalues found. \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
          The real parts of each scalar alpha defining an eigenvalue
          of GNEP. \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
          The imaginary parts of each scalar alpha defining an
          eigenvalue of GNEP. \n
          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
          positive, then the j-th and (j+1)-st eigenvalues are a
          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j). \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          The scalars beta that define the eigenvalues of GNEP. \n
          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
          beta = BETA(j) represent the j-th eigenvalue of the matrix
          pair (A,B), in one of the forms lambda = alpha/beta or
          mu = beta/alpha.  Since either lambda or mu may overflow,
          they should not, in general, be computed. \n
 * @param[in, out] QC
          QC is REAL array, dimension (LDQC, NW) \n
 * @param[in] LDQC
          LDQC is INTEGER \n
 * @param[in, out] ZC
          ZC is REAL array, dimension (LDZC, NW) \n
 * @param[in] LDZ
          LDZ is INTEGER \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[in] REC
          REC is INTEGER \n
          REC indicates the current recursion level. Should be set
          to 0 on first call. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void laqz3(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi,
               integer *nw, T *a, integer *lda, T *b, integer *ldb, T *q, integer *ldq, T *z,
               integer *ldz, integer *ns, integer *nd, T *alphar, T *alphai, T *beta, T *qc,
               integer *ldqc, T *zc, integer *ldzc, T *work, integer *lwork, integer *rec,
               integer *info)
    {
        laqz3(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, ns, nd, alphar,
              alphai, beta, qc, ldqc, zc, ldzc, work, lwork, rec, info);
    }
    template <typename T>
    void laqz3(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi,
               integer *nshifts, integer *nblock_desired, T *alpha, T *beta, T *a, integer *lda,
               T *b, integer *ldb, T *q, integer *ldq, T *z, integer *ldz, T *qc, integer *ldqc,
               T *zc, integer *ldzc, T *work, integer *lwork, integer *info)
    {
        laqz3(ilschur, ilq, ilz, n, ilo, ihi, nshifts, nblock_desired, alpha, beta, a, lda, b, ldb,
              q, ldq, z, ldz, qc, ldqc, zc, ldzc, work, lwork, info);
    }
    /** @}*/ // end of laqz3

    /** @defgroup laqz4 laqz4
     * @ingroup NYS
     * @{
     */
    /*! @brief LAQZ4 Executes a single multishift QZ sweep

 * @details
 * \b Purpose:
    \verbatim
    LAQZ4 Executes a single multishift QZ sweep
    \endverbatim

 * @param[in] ILSCHUR
          ILSCHUR is LOGICAL \n
          Determines whether or not to update the full Schur form \n
 * @param[in] ILQ
          ILQ is LOGICAL \n
          Determines whether or not to update the matrix Q \n
 * @param[in] ILZ
          ILZ is LOGICAL \n
          Determines whether or not to update the matrix Z \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A, B, Q, and Z.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
 * @param[in] NSHIFTS
          NSHIFTS is INTEGER \n
          The desired number of shifts to use \n
 * @param[in] NBLOCK_DESIRED
          NBLOCK_DESIRED is INTEGER \n
          The desired size of the computational windows \n
 * @param[in] SR
          SR is REAL array. SR contains
          the real parts of the shifts to use. \n
 * @param[in] SI
          SI is REAL array. SI contains
          the imaginary parts of the shifts to use. \n
 * @param[in] SS
          SS is REAL array. SS contains
          the scale of the shifts to use. \n
 * @param[in, out] A
          A is REAL array, dimension (LDA, N) \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max( 1, N ). \n
 * @param[in, out] B
          B is REAL array, dimension (LDB,N) \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max( 1, N ). \n
 * @param[in, out] Q
          Q is REAL array, dimension (LDQ,N) \n
 * @param[in] LDQ
          LDQ is INTEGER \n
 * @param[in, out] Z
          Z is REAL array, dimension (LDZ,N) \n
 * @param[in] LDZ
          LDZ is INTEGER \n
 * @param[in, out] QC
          QC is REAL array, dimension (LDQC, NBLOCK_DESIRED) \n
 * @param[in] LDQC
          LDQC is INTEGER \n
 * @param[in, out] ZC
          ZC is REAL array, dimension (LDZC, NBLOCK_DESIRED) \n
 * @param[in] LDZC
          LDZC is INTEGER \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void laqz4(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi,
               integer *nshifts, integer *nblock_desired, T *sr, T *si, T *ss, T *a, integer *lda,
               T *b, integer *ldb, T *q, integer *ldq, T *z, integer *ldz, T *qc, integer *ldqc,
               T *zc, integer *ldzc, T *work, integer *lwork, integer *info)
    {
        laqz4(ilschur, ilq, ilz, n, ilo, ihi, nshifts, nblock_desired, sr, si, ss, a, lda, b, ldb,
              q, ldq, z, ldz, qc, ldqc, zc, ldzc, work, lwork, info);
    }
    /** @}*/ // end of laqz4

    /** @defgroup ggbal ggbal
     * @ingroup NYS
     * @{
     */
    /*! @brief GGBAL balances a pair of general real matrices (A,B)
 * @details
 * \b Purpose:
    \verbatim
     SGGBAL balances a pair of general real matrices (A,B).  This
     involves, first, permuting A and B by similarity transformations to
     isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N
     elements on the diagonal; and second, applying a diagonal similarity
     transformation to rows and columns ILO to IHI to make the rows
     and columns as close in norm as possible. Both steps are optional.

     Balancing may reduce the 1-norm of the matrices, and improve the
     accuracy of the computed eigenvalues and/or eigenvectors in the
     generalized eigenvalue problem A*x = lambda*B*x.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          Specifies the operations to be performed on A and B: \n
          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
                  and RSCALE(I) = 1.0 for i = 1,...,N. \n
          = 'P':  permute only; \n
          = 'S':  scale only; \n
          = 'B':  both permute and scale. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the input matrix A. \n
          On exit,  A is overwritten by the balanced matrix.
          If JOB = 'N', A is not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the input matrix B. \n
          On exit,  B is overwritten by the balanced matrix.
          If JOB = 'N', B is not referenced. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,N). \n
 * @param[out] ILO
          ILO is INTEGER \n
 * @param[out] IHI
          IHI is INTEGER \n
          ILO and IHI are set to integers such that on exit
          A(i,j) = 0 and B(i,j) = 0 if i > j and
          j = 1,...,ILO-1 or i = IHI+1,...,N.
          If JOB = 'N' or 'S', ILO = 1 and IHI = N. \n
 * @param[out] LSCALE
          LSCALE is REAL array, dimension (N) \n
          Details of the permutations and scaling factors applied
          to the left side of A and B.  If P(j) is the index of the
          row interchanged with row j, and D(j)
          is the scaling factor applied to row j, then \n
            LSCALE(j) = P(j)    for J = 1,...,ILO-1 \n
                      = D(j)    for J = ILO,...,IHI \n
                      = P(j)    for J = IHI+1,...,N. \n
          The order in which the interchanges are made is N to IHI+1,
          then 1 to ILO-1. \n
 * @param[out] RSCALE
          RSCALE is REAL array, dimension (N) \n
          Details of the permutations and scaling factors applied
          to the right side of A and B.  If P(j) is the index of the
          column interchanged with column j, and D(j)
          is the scaling factor applied to column j, then \n
            LSCALE(j) = P(j)    for J = 1,...,ILO-1 \n
                      = D(j)    for J = ILO,...,IHI \n
                      = P(j)    for J = IHI+1,...,N. \n
          The order in which the interchanges are made is N to IHI+1,
          then 1 to ILO-1. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (lwork) \n
          lwork must be at least fla_max(1,6*N) when JOB = 'S' or 'B', and
          at least 1 when JOB = 'N' or 'P'. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void ggbal(char *job, integer *n, T *a, integer *lda, T *b, integer *ldb, integer *ilo,
               integer *ihi, T *lscale, T *rscale, T *work, integer *info)
    {
        ggbal(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
    }
    template <typename T, typename Ta>
    void ggbal(char *job, integer *n, T *a, integer *lda, T *b, integer *ldb, integer *ilo,
               integer *ihi, Ta *lscale, Ta *rscale, Ta *work, integer *info)
    {
        ggbal(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
    }
    /** @}*/ // end of ggbal

    /** @defgroup gghrd gghrd
     * @ingroup NYS
     * @{
     */
    /*! @brief GGHRD reduces a pair of real matrices (A,B) to generalized upper Hessenberg \n
    form using orthogonal transformations
 * @details
 * \b Purpose:
    \verbatim
     GGHRD reduces a pair of real matrices (A,B) to generalized upper
     Hessenberg form using orthogonal transformations, where A is a
     general matrix and B is upper triangular.  The form of the
     generalized eigenvalue problem is
        A*x = lambda*B*x,
     and B is typically made upper triangular by computing its QR
     factorization and moving the orthogonal matrix Q to the left side
     of the equation.

     This subroutine simultaneously reduces A to a Hessenberg matrix H:
        Q**T*A*Z = H
     and transforms B to another upper triangular matrix T:
        Q**T*B*Z = T
     in order to reduce the problem to its standard form
        H*y = lambda*T*y
     where y = Z**T*x.

     The orthogonal matrices Q and Z are determined as products of Givens
     rotations.  They may either be formed explicitly, or they may be
     postmultiplied into input matrices Q1 and Z1, so that
          Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T
          Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T

     If Q1 is the orthogonal matrix from the QR factorization of B in the
     original equation A*x = lambda*B*x, then SGGHD3 reduces the original
     problem to generalized Hessenberg form.

     This is a blocked variant of SGGHRD, using matrix-matrix
     multiplications for parts of the computation to enhance performance.
    \endverbatim

 * @param[in] COMPQ
          COMPQ is CHARACTER*1 \n
          = 'N': do not compute Q; \n
          = 'I': Q is initialized to the unit matrix, and the
                 orthogonal matrix Q is returned; \n
          = 'V': Q must contain an orthogonal matrix Q1 on entry,
                 and the product Q1*Q is returned. \n
 * @param[in] COMPZ
          COMPZ is CHARACTER*1 \n
          = 'N': do not compute Z; \n
          = 'I': Z is initialized to the unit matrix, and the
                 orthogonal matrix Z is returned; \n
          = 'V': Z must contain an orthogonal matrix Z1 on entry,
                 and the product Z1*Z is returned. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
 \n
          ILO and IHI mark the rows and columns of A which are to be
          reduced.  It is assumed that A is already upper triangular
          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
          normally set by a previous call to SGGBAL; otherwise they
          should be set to 1 and N respectively. \n
          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the N-by-N general matrix to be reduced. \n
          On exit, the upper triangle and the first subdiagonal of A
          are overwritten with the upper Hessenberg matrix H, and the
          rest is set to zero. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the N-by-N upper triangular matrix B. \n
          On exit, the upper triangular matrix T = Q**T B Z.  The
          elements below the diagonal are set to zero. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ, N) \n
          On entry, if COMPQ = 'V', the orthogonal matrix Q1,
          typically from the QR factorization of B. \n
          On exit, if COMPQ='I', the orthogonal matrix Q, and if
          COMPQ = 'V', the product Q1*Q.
          Not referenced if COMPQ='N'. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. \n
          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ, N) \n
          On entry, if COMPZ = 'V', the orthogonal matrix Z1. \n
          On exit, if COMPZ='I', the orthogonal matrix Z, and if
          COMPZ = 'V', the product Z1*Z.
          Not referenced if COMPZ='N'. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z. \n
          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void gghrd(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, T *a, integer *lda,
               T *b, integer *ldb, T *q, integer *ldq, T *z, integer *ldz, integer *info)
    {
        gghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
    }
    /** @}*/ // end of gghrd

    /** @defgroup gghd3 gghd3
     * @ingroup NYS
     * @{
     */
    /*! @brief GGHD3 reduces a pair of real matrices (A,B) to generalized upper Hessenberg \n
    form using orthogonal transformations

 * @details
 * \b Purpose:
    \verbatim
     GGHD3 reduces a pair of real matrices (A,B) to generalized upper
     Hessenberg form using orthogonal transformations, where A is a
     general matrix and B is upper triangular.  The form of the
     generalized eigenvalue problem is
        A*x = lambda*B*x,
     and B is typically made upper triangular by computing its QR
     factorization and moving the orthogonal matrix Q to the left side
     of the equation.

     This subroutine simultaneously reduces A to a Hessenberg matrix H:
        Q**T*A*Z = H
     and transforms B to another upper triangular matrix T:
        Q**T*B*Z = T
     in order to reduce the problem to its standard form
        H*y = lambda*T*y
     where y = Z**T*x.

     The orthogonal matrices Q and Z are determined as products of Givens
     rotations.  They may either be formed explicitly, or they may be
     postmultiplied into input matrices Q1 and Z1, so that
          Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T
          Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T

     If Q1 is the orthogonal matrix from the QR factorization of B in the
     original equation A*x = lambda*B*x, then SGGHD3 reduces the original
     problem to generalized Hessenberg form.

     This is a blocked variant of SGGHRD, using matrix-matrix
     multiplications for parts of the computation to enhance performance.
    \endverbatim

 * @param[in] COMPQ
          COMPQ is CHARACTER*1 \n
          = 'N': do not compute Q; \n
          = 'I': Q is initialized to the unit matrix, and the
                 orthogonal matrix Q is returned; \n
          = 'V': Q must contain an orthogonal matrix Q1 on entry,
                 and the product Q1*Q is returned. \n
 * @param[in] COMPZ
          COMPZ is CHARACTER*1 \n
          = 'N': do not compute Z; \n
          = 'I': Z is initialized to the unit matrix, and the
                 orthogonal matrix Z is returned; \n
          = 'V': Z must contain an orthogonal matrix Z1 on entry,
                 and the product Z1*Z is returned. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER
 \n
          ILO and IHI mark the rows and columns of A which are to be
          reduced.  It is assumed that A is already upper triangular
          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
          normally set by a previous call to SGGBAL; otherwise they
          should be set to 1 and N respectively. \n
          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the N-by-N general matrix to be reduced. \n
          On exit, the upper triangle and the first subdiagonal of A
          are overwritten with the upper Hessenberg matrix H, and the
          rest is set to zero. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the N-by-N upper triangular matrix B. \n
          On exit, the upper triangular matrix T = Q**T B Z.  The
          elements below the diagonal are set to zero. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ, N) \n
          On entry, if COMPQ = 'V', the orthogonal matrix Q1,
          typically from the QR factorization of B. \n
          On exit, if COMPQ='I', the orthogonal matrix Q, and if
          COMPQ = 'V', the product Q1*Q. \n
          Not referenced if COMPQ='N'. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. \n
          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ, N) \n
          On entry, if COMPZ = 'V', the orthogonal matrix Z1. \n
          On exit, if COMPZ='I', the orthogonal matrix Z, and if
          COMPZ = 'V', the product Z1*Z. \n
          Not referenced if COMPZ='N'. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z. \n
          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (LWORK) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of the array WORK.  LWORK >= 1. \n
          For optimum performance LWORK >= 6*N*NB, where NB is the
          optimal blocksize. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void gghd3(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, T *a, integer *lda,
               T *b, integer *ldb, T *q, integer *ldq, T *z, integer *ldz, T *work, integer *lwork,
               integer *info)
    {
        gghd3(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
    }
    /** @}*/ // end of gghd3

    /** @defgroup hgeqz hgeqz
     * @ingroup NYS
     * @{
     */
    /*! @brief HGEQZ computes the eigenvalues of a real matrix pair (H,T).

 * @details
 * \b Purpose:
    \verbatim
     HGEQZ computes the eigenvalues of a real matrix pair (H,T),
     where H is an upper Hessenberg matrix and T is upper triangular,
     using the double-shift QZ method.
     Matrix pairs of this type are produced by the reduction to
     generalized upper Hessenberg form of a real matrix pair (A,B):

        A = Q1*H*Z1**T,  B = Q1*T*Z1**T,

     as computed by SGGHRD.

     If JOB='S', then the Hessenberg-triangular pair (H,T) is
     also reduced to generalized Schur form,

        H = Q*S*Z**T,  T = Q*P*Z**T,

     where Q and Z are orthogonal matrices, P is an upper triangular
     matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2
     diagonal blocks.

     The 1-by-1 blocks correspond to real eigenvalues of the matrix pair
     (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of
     eigenvalues.

     Additionally, the 2-by-2 upper triangular diagonal blocks of P
     corresponding to 2-by-2 blocks of S are reduced to positive diagonal
     form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,
     P(j,j) > 0, and P(j+1,j+1) > 0.

     Optionally, the orthogonal matrix Q from the generalized Schur
     factorization may be postmultiplied into an input matrix Q1, and the
     orthogonal matrix Z may be postmultiplied into an input matrix Z1.
     If Q1 and Z1 are the orthogonal matrices from SGGHRD that reduced
     the matrix pair (A,B) to generalized upper Hessenberg form, then the
     output matrices Q1*Q and Z1*Z are the orthogonal factors from the
     generalized Schur factorization of (A,B):

        A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.

     To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,
     of (A,B)) are computed as a pair of values (alpha,beta), where alpha is
     complex and beta real.
     If beta is nonzero, lambda = alpha / beta is an eigenvalue of the
     generalized nonsymmetric eigenvalue problem (GNEP)
        A*x = lambda*B*x
     and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
     alternate form of the GNEP
        mu*A*y = B*y.
     Real eigenvalues can be read directly from the generalized Schur
     form:
       alpha = S(i,i), beta = P(i,i).

     Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
          Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
          pp. 241--256.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          = 'E': Compute eigenvalues only; \n
          = 'S': Compute eigenvalues and the Schur form. \n
 * @param[in] COMPQ
          COMPQ is CHARACTER*1 \n
          = 'N': Left Schur vectors (Q) are not computed; \n
          = 'I': Q is initialized to the unit matrix and the matrix Q
                 of left Schur vectors of (H,T) is returned; \n
          = 'V': Q must contain an orthogonal matrix Q1 on entry and
                 the product Q1*Q is returned. \n
 * @param[in] COMPZ
          COMPZ is CHARACTER*1 \n
          = 'N': Right Schur vectors (Z) are not computed; \n
          = 'I': Z is initialized to the unit matrix and the matrix Z
                 of right Schur vectors of (H,T) is returned; \n
          = 'V': Z must contain an orthogonal matrix Z1 on entry and
                 the product Z1*Z is returned. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices H, T, Q, and Z.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          ILO and IHI mark the rows and columns of H which are in
          Hessenberg form.  It is assumed that A is already upper
          triangular in rows and columns 1:ILO-1 and IHI+1:N.
          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0. \n
 * @param[in,out] H
          H is REAL array, dimension (LDH, N) \n
          On entry, the N-by-N upper Hessenberg matrix H. \n
          On exit, if JOB = 'S', H contains the upper quasi-triangular
          matrix S from the generalized Schur factorization.
          If JOB = 'E', the diagonal blocks of H match those of S, but
          the rest of H is unspecified. \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the array H.  LDH >= fla_max( 1, N). \n
 * @param[in,out] T
          T is REAL array, dimension (LDT, N) \n
          On entry, the N-by-N upper triangular matrix T. \n
          On exit, if JOB = 'S', T contains the upper triangular
          matrix P from the generalized Schur factorization;
          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S
          are reduced to positive diagonal form, i.e., if H(j+1,j) is
          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and
          T(j+1,j+1) > 0. \n
          If JOB = 'E', the diagonal blocks of T match those of P, but
          the rest of T is unspecified. \n
 * @param[in] LDT
          LDT is INTEGER \n
          The leading dimension of the array T.  LDT >= fla_max( 1, N). \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
          The real parts of each scalar alpha defining an eigenvalue
          of GNEP. \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
          The imaginary parts of each scalar alpha defining an
          eigenvalue of GNEP. \n
          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
          positive, then the j-th and (j+1)-st eigenvalues are a
          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j). \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          The scalars beta that define the eigenvalues of GNEP.
          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
          beta = BETA(j) represent the j-th eigenvalue of the matrix
          pair (A,B), in one of the forms lambda = alpha/beta or
          mu = beta/alpha.  Since either lambda or mu may overflow,
          they should not, in general, be computed. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ, N) \n
          On entry, if COMPQ = 'V', the orthogonal matrix Q1 used in
          the reduction of (A,B) to generalized Hessenberg form.
          On exit, if COMPQ = 'I', the orthogonal matrix of left Schur
          vectors of (H,T), and if COMPQ = 'V', the orthogonal matrix
          of left Schur vectors of (A,B).
          Not referenced if COMPQ = 'N'. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= 1.
          If COMPQ='V' or 'I', then LDQ >= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ, N) \n
          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in
          the reduction of (A,B) to generalized Hessenberg form.
          On exit, if COMPZ = 'I', the orthogonal matrix of
          right Schur vectors of (H,T), and if COMPZ = 'V', the
          orthogonal matrix of right Schur vectors of (A,B).
          Not referenced if COMPZ = 'N'. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1.
          If COMPZ='V' or 'I', then LDZ >= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER
          The dimension of the array WORK.  LWORK >= fla_max(1,N). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
                     in Schur form, but ALPHAR(i), ALPHAI(i), and
                     BETA(i), i=INFO+1,...,N should be correct. \n
          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
                     in Schur form, but ALPHAR(i), ALPHAI(i), and
                     BETA(i), i=INFO-N+1,...,N should be correct. \n

 *  * */
    template <typename T>
    void hgeqz(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *t, integer *ldt, T *alphar, T *alphai, T *beta, T *q, integer *ldq,
               T *z, integer *ldz, T *work, integer *lwork, integer *info)
    {
        hgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz,
              work, lwork, info);
    }
    template <typename T, typename Ta>
    void hgeqz(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi, T *h,
               integer *ldh, T *t, integer *ldt, T *alpha, T *beta, T *q, integer *ldq, T *z,
               integer *ldz, T *work, integer *lwork, Ta *rwork, integer *info)
    {
        hgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz, work,
              lwork, rwork, info);
    }
    /** @}*/ // end of hgeqz

    /** @defgroup ggbak ggbak
     * @ingroup NYS
     * @{
     */
    /*! @brief GGBAK forms the right or left eigenvectors of a real generalized eigenvalue problem
 * @details
 * \b Purpose:
    \verbatim
     GGBAK forms the right or left eigenvectors of a real generalized
     eigenvalue problem A*x = lambda*B*x, by backward transformation on
     the computed eigenvectors of the balanced pair of matrices output by
     GGBAL.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          Specifies the type of backward transformation required: \n
          = 'N':  do nothing, return immediately; \n
          = 'P':  do backward transformation for permutation only; \n
          = 'S':  do backward transformation for scaling only; \n
          = 'B':  do backward transformations for both permutation and
                  scaling. \n
          JOB must be the same as the argument JOB supplied to SGGBAL. \n
 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'R':  V contains right eigenvectors; \n
          = 'L':  V contains left eigenvectors. \n
 * @param[in] N
          N is INTEGER \n
          The number of rows of the matrix V.  N >= 0. \n
 * @param[in] ILO
          ILO is INTEGER \n
 * @param[in] IHI
          IHI is INTEGER \n
          The integers ILO and IHI determined by SGGBAL. \n
          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. \n
 * @param[in] LSCALE
          LSCALE is REAL array, dimension (N) \n
          Details of the permutations and/or scaling factors applied
          to the left side of A and B, as returned by SGGBAL. \n
 * @param[in] RSCALE
          RSCALE is REAL array, dimension (N) \n
          Details of the permutations and/or scaling factors applied
          to the right side of A and B, as returned by SGGBAL. \n
 * @param[in] M
          M is INTEGER \n
          The number of columns of the matrix V.  M >= 0. \n
 * @param[in,out] V
          V is REAL array, dimension (LDV,M) \n
          On entry, the matrix of right or left eigenvectors to be
          transformed, as returned by STGEVC. \n
          On exit, V is overwritten by the transformed eigenvectors. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the matrix V. LDV >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void ggbak(char *job, char *side, integer *n, integer *ilo, integer *ihi, T *lscale, T *rscale,
               integer *m, T *v, integer *ldv, integer *info)
    {
        ggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info);
    }
    template <typename T, typename Ta>
    void ggbak(char *job, char *side, integer *n, integer *ilo, integer *ihi, Ta *lscale,
               Ta *rscale, integer *m, T *v, integer *ldv, integer *info)
    {
        ggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info);
    }
    /** @}*/ // end of ggbak

    /** @defgroup tgsen tgsen
     * @ingroup NYS
     * @{
     */
    /*! @brief TGSEN reorders the generalized real Schur decomposition of a real matrix pair

 * @details
 * \b Purpose:
    \verbatim
     TGSEN reorders the generalized real Schur decomposition of a real
     matrix pair (A, B) (in terms of an orthonormal equivalence trans-
     formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues
     appears in the leading diagonal blocks of the upper quasi-triangular
     matrix A and the upper triangular B. The leading columns of Q and
     Z form orthonormal bases of the corresponding left and right eigen-
     spaces (deflating subspaces). (A, B) must be in generalized real
     Schur canonical form (as returned by SGGES), i.e. A is block upper
     triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper
     triangular.

     STGSEN also computes the generalized eigenvalues

                 w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)

     of the reordered matrix pair (A, B).

     Optionally, STGSEN computes the estimates of reciprocal condition
     numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     the selected cluster and the eigenvalues outside the cluster, resp.,
     and norms of "projections" onto left and right eigenspaces w.r.t.
     the selected cluster in the (1,1)-block.
    \endverbatim

 * @param[in] IJOB
          IJOB is INTEGER \n
          Specifies whether condition numbers are required for the
          cluster of eigenvalues (PL and PR) or the deflating subspaces
          (Difu and Difl): \n
           =0: Only reorder w.r.t. SELECT. No extras. \n
           =1: Reciprocal of norms of "projections" onto left and right
               eigenspaces w.r.t. the selected cluster (PL and PR). \n
           =2: Upper bounds on Difu and Difl. F-norm-based estimate
               (DIF(1:2)). \n
           =3: Estimate of Difu and Difl. 1-norm-based estimate
               (DIF(1:2)).
               About 5 times as expensive as IJOB = 2. \n
           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic
               version to get it all. \n
           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above) \n
 * @param[in] WANTQ
          WANTQ is LOGICAL \n
          .TRUE. : update the left transformation matrix Q; \n
          .FALSE.: do not update Q. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          .TRUE. : update the right transformation matrix Z; \n
          .FALSE.: do not update Z. \n
 * @param[in] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          SELECT specifies the eigenvalues in the selected cluster.
          To select a real eigenvalue w(j), SELECT(j) must be set to
          .TRUE.. To select a complex conjugate pair of eigenvalues
          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
          either SELECT(j) or SELECT(j+1) or both must be set to
          .TRUE.; a complex conjugate pair of eigenvalues must be
          either both included in the cluster or both excluded. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension(LDA,N) \n
          On entry, the upper quasi-triangular matrix A, with (A, B) in
          generalized real Schur canonical form.
          On exit, A is overwritten by the reordered matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension(LDB,N) \n
          On entry, the upper triangular matrix B, with (A, B) in
          generalized real Schur canonical form. \n
          On exit, B is overwritten by the reordered matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,N). \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (N) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i
          and BETA(j),j=1,...,N  are the diagonals of the complex Schur
          form (S,T) that would result if the 2-by-2 diagonal blocks of
          the real generalized Schur form of (A,B) were further reduced
          to triangular form using complex unitary transformations.
          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
          positive, then the j-th and (j+1)-st eigenvalues are a
          complex conjugate pair, with ALPHAI(j+1) negative. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix. \n
          On exit, Q has been postmultiplied by the left orthogonal
          transformation matrix which reorder (A, B); The leading M
          columns of Q form orthonormal bases for the specified pair of
          left eigenspaces (deflating subspaces). \n
          If WANTQ = .FALSE., Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= 1; \n
          and if WANTQ = .TRUE., LDQ >= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,N) \n
          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix. \n
          On exit, Z has been postmultiplied by the left orthogonal
          transformation matrix which reorder (A, B); The leading M
          columns of Z form orthonormal bases for the specified pair of
          left eigenspaces (deflating subspaces). \n
          If WANTZ = .FALSE., Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z. LDZ >= 1;
          If WANTZ = .TRUE., LDZ >= N. \n
 * @param[out] M
          M is INTEGER \n
          The dimension of the specified pair of left and right eigen-
          spaces (deflating subspaces). 0 <= M <= N. \n
 * @param[out] PL
          PL is REAL \n
 * @param[out] PR
          PR is REAL \n
          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the
          reciprocal of the norm of "projections" onto left and right
          eigenspaces with respect to the selected cluster. \n
          0 < PL, PR <= 1. \n
          If M = 0 or M = N, PL = PR  = 1. \n
          If IJOB = 0, 2 or 3, PL and PR are not referenced. \n
 * @param[out] DIF
          DIF is REAL array, dimension (2). \n
          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
          estimates of Difu and Difl. \n
          If M = 0 or N, DIF(1:2) = F-norm([A, B]). \n
          If IJOB = 0 or 1, DIF is not referenced. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >=  4*N+16. \n
          If IJOB = 1, 2 or 4, LWORK >= MAX(4*N+16, 2*M*(N-M)). \n
          If IJOB = 3 or 5, LWORK >= MAX(4*N+16, 4*M*(N-M)). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. LIWORK >= 1. \n
          If IJOB = 1, 2 or 4, LIWORK >=  N+6. \n
          If IJOB = 3 or 5, LIWORK >= MAX(2*M*(N-M), N+6). \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the IWORK array,
          returns this value as the first entry of the IWORK array, and
          no error message related to LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
            =0: Successful exit. \n
            <0: If INFO = -i, the i-th argument had an illegal value. \n
            =1: Reordering of (A, B) failed because the transformed
                matrix pair (A, B) would be too far from generalized
                Schur form; the problem is very ill-conditioned.
                (A, B) may have been partially reordered.
                If requested, 0 is returned in DIF(*), PL and PR. \n

 *  * */
    template <typename T>
    void tgsen(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n, T *a,
               integer *lda, T *b, integer *ldb, T *alphar, T *alphai, T *beta, T *q, integer *ldq,
               T *z, integer *ldz, integer *m, T *pl, T *pr, T *dif, T *work, integer *lwork,
               integer *iwork, integer *liwork, integer *info)
    {
        tgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz,
              m, pl, pr, dif, work, lwork, iwork, liwork, info);
    }
    template <typename T, typename Ta>
    void tgsen(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n, T *a,
               integer *lda, T *b, integer *ldb, T *alpha, T *beta, T *q, integer *ldq, T *z,
               integer *ldz, integer *m, Ta *pl, Ta *pr, Ta *dif, T *work, integer *lwork,
               integer *iwork, integer *liwork, integer *info)
    {
        tgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, m, pl, pr,
              dif, work, lwork, iwork, liwork, info);
    }
    /** @}*/ // end of tgsen

    /** @defgroup tgsna tgsna
     * @ingroup NYS
     * @{
     */
    /*! @brief TGSNA estimates reciprocal condition numbers for specified  \n
     eigenvalues and/or eigenvectors of a matrix pair
 * @details
 * \b Purpose:
    \verbatim
     TGSNA estimates reciprocal condition numbers for specified
     eigenvalues and/or eigenvectors of a matrix pair (A, B) in
     generalized real Schur canonical form (or of any matrix pair
     (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where
     Z**T denotes the transpose of Z.

     (A, B) must be in generalized real Schur form (as returned by SGGES),
     i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal
     blocks. B is upper triangular.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          Specifies whether condition numbers are required for
          eigenvalues (S) or eigenvectors (DIF): \n
          = 'E': for eigenvalues only (S); \n
          = 'V': for eigenvectors only (DIF); \n
          = 'B': for both eigenvalues and eigenvectors (S and DIF). \n
 * @param[in] HOWMNY
          HOWMNY is CHARACTER*1 \n
          = 'A': compute condition numbers for all eigenpairs; \n
          = 'S': compute condition numbers for selected eigenpairs
                 specified by the array SELECT. \n
 * @param[in] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
          condition numbers are required. To select condition numbers
          for the eigenpair corresponding to a real eigenvalue w(j),
          SELECT(j) must be set to .TRUE.. To select condition numbers
          corresponding to a complex conjugate pair of eigenvalues w(j)
          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
          set to .TRUE.. \n
          If HOWMNY = 'A', SELECT is not referenced. \n
 * @param[in] N
          N is INTEGER \n
          The order of the square matrix pair (A, B). N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The upper quasi-triangular matrix A in the pair (A,B). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[in] B
          B is REAL array, dimension (LDB,N) \n
          The upper triangular matrix B in the pair (A,B). \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,N). \n
 * @param[in] VL
          VL is REAL array, dimension (LDVL,M) \n
          If JOB = 'E' or 'B', VL must contain left eigenvectors of
          (A, B), corresponding to the eigenpairs specified by HOWMNY
          and SELECT. The eigenvectors must be stored in consecutive
          columns of VL, as returned by STGEVC. \n
          If JOB = 'V', VL is not referenced. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of the array VL. LDVL >= 1.
          If JOB = 'E' or 'B', LDVL >= N. \n
 * @param[in] VR
          VR is REAL array, dimension (LDVR,M) \n
          If JOB = 'E' or 'B', VR must contain right eigenvectors of
          (A, B), corresponding to the eigenpairs specified by HOWMNY
          and SELECT. The eigenvectors must be stored in consecutive
          columns ov VR, as returned by STGEVC. \n
          If JOB = 'V', VR is not referenced. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR. LDVR >= 1.
          If JOB = 'E' or 'B', LDVR >= N. \n
 * @param[out] S
          S is REAL array, dimension (MM) \n
          If JOB = 'E' or 'B', the reciprocal condition numbers of the
          selected eigenvalues, stored in consecutive elements of the
          array. For a complex conjugate pair of eigenvalues two
          consecutive elements of S are set to the same value. Thus
          S(j), DIF(j), and the j-th columns of VL and VR all
          correspond to the same eigenpair (but not in general the
          j-th eigenpair, unless all eigenpairs are selected). \n
          If JOB = 'V', S is not referenced. \n
 * @param[out] DIF
          DIF is REAL array, dimension (MM) \n
          If JOB = 'V' or 'B', the estimated reciprocal condition
          numbers of the selected eigenvectors, stored in consecutive
          elements of the array. For a complex eigenvector two
          consecutive elements of DIF are set to the same value. If
          the eigenvalues cannot be reordered to compute DIF(j), DIF(j)
          is set to 0; this can only occur when the true value would be
          very small anyway. \n
          If JOB = 'E', DIF is not referenced. \n
 * @param[in] MM
          MM is INTEGER \n
          The number of elements in the arrays S and DIF. MM >= M. \n
 * @param[out] M
          M is INTEGER \n
          The number of elements of the arrays S and DIF used to store
          the specified condition numbers; for each selected real
          eigenvalue one element is used, and for each selected complex
          conjugate pair of eigenvalues, two elements are used.
          If HOWMNY = 'A', M is set to N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,N).
          If JOB = 'V' or 'B' LWORK >= 2*N*(N+2)+16. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N + 6)
          If JOB = 'E', IWORK is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          =0: Successful exit \n
          <0: If INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tgsna(char *job, char *howmny, logical *select, integer *n, T *a, integer *lda, T *b,
               integer *ldb, T *vl, integer *ldvl, T *vr, integer *ldvr, T *s, T *dif, integer *mm,
               integer *m, T *work, integer *lwork, integer *iwork, integer *info)
    {
        tgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work,
              lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void tgsna(char *job, char *howmny, logical *select, integer *n, T *a, integer *lda, T *b,
               integer *ldb, T *vl, integer *ldvl, T *vr, integer *ldvr, Ta *s, Ta *dif,
               integer *mm, integer *m, T *work, integer *lwork, integer *iwork, integer *info)
    {
        tgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work,
              lwork, iwork, info);
    }
    /** @}*/ // end of tgsna

    /** @defgroup tgsyl tgsyl
     * @ingroup NYS
     * @{
     */
    /*! @brief TGSYL solves the generalized Sylvester equation

 * @details
 * \b Purpose:
    \verbatim
     TGSYL solves the generalized Sylvester equation:

                 A * R - L * B = scale * C                 (1)
                 D * R - L * E = scale * F

     where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     respectively, with real entries. (A, D) and (B, E) must be in
     generalized (real) Schur canonical form, i.e. A, B are upper quasi
     triangular and D, E are upper triangular.

     The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     scaling factor chosen to avoid overflow.

     In matrix notation (1) is equivalent to solve  Zx = scale b, where
     Z is defined as

                Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)
                    [ kron(In, D)  -kron(E**T, Im) ].

     Here Ik is the identity matrix of size k and X**T is the transpose of
     X. kron(X, Y) is the Kronecker product between the matrices X and Y.

     If TRANS = 'T', STGSYL solves the transposed system Z**T*y = scale*b,
     which is equivalent to solve for R and L in

                 A**T * R + D**T * L = scale * C           (3)
                 R * B**T + L * E**T = scale * -F

     This case (TRANS = 'T') is used to compute an one-norm-based estimate
     of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     and (B,E), using SLACON.

     If IJOB >= 1, STGSYL computes a Frobenius norm-based estimate
     of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     reciprocal of the smallest singular value of Z. See [1-2] for more
     information.

     This is a level 3 BLAS algorithm.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': solve the generalized Sylvester equation (1).
          = 'T': solve the 'transposed' system (3). \n
 * @param[in] IJOB
          IJOB is INTEGER \n
          Specifies what kind of functionality to be performed. \n
          = 0: solve (1) only. \n
          = 1: The functionality of 0 and 3. \n
          = 2: The functionality of 0 and 4. \n
          = 3: Only an estimate of Dif[(A,D), (B,E)] is computed.
               (look ahead strategy IJOB  = 1 is used). \n
          = 4: Only an estimate of Dif[(A,D), (B,E)] is computed.
               ( SGECON on sub-systems is used). \n
          Not referenced if TRANS = 'T'. \n
 * @param[in] M
          M is INTEGER \n
          The order of the matrices A and D, and the row dimension of
          the matrices C, F, R and L. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices B and E, and the column dimension
          of the matrices C, F, R and L. \n
 * @param[in] A
          A is REAL array, dimension (LDA, M) \n
          The upper quasi triangular matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1, M). \n
 * @param[in] B
          B is REAL array, dimension (LDB, N) \n
          The upper quasi triangular matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1, N). \n
 * @param[in,out] C
          C is REAL array, dimension (LDC, N) \n
          On entry, C contains the right-hand-side of the first matrix
          equation in (1) or (3). \n
          On exit, if IJOB = 0, 1 or 2, C has been overwritten by
          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,
          the solution achieved during the computation of the
          Dif-estimate. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1, M). \n
 * @param[in] D
          D is REAL array, dimension (LDD, M) \n
          The upper triangular matrix D. \n
 * @param[in] LDD
          LDD is INTEGER \n
          The leading dimension of the array D. LDD >= fla_max(1, M). \n
 * @param[in] E
          E is REAL array, dimension (LDE, N) \n
          The upper triangular matrix E. \n
 * @param[in] LDE
          LDE is INTEGER \n
          The leading dimension of the array E. LDE >= fla_max(1, N). \n
 * @param[in,out] F
          F is REAL array, dimension (LDF, N) \n
          On entry, F contains the right-hand-side of the second matrix
          equation in (1) or (3). \n
          On exit, if IJOB = 0, 1 or 2, F has been overwritten by
          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,
          the solution achieved during the computation of the
          Dif-estimate. \n
 * @param[in] LDF
          LDF is INTEGER \n
          The leading dimension of the array F. LDF >= fla_max(1, M). \n
 * @param[out] DIF
          DIF is REAL \n
          On exit DIF is the reciprocal of a lower bound of the
          reciprocal of the Dif-function, i.e. DIF is an upper bound of
          Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2).
          IF IJOB = 0 or TRANS = 'T', DIF is not touched. \n
 * @param[out] SCALE
          SCALE is REAL \n
          On exit SCALE is the scaling factor in (1) or (3).
          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,
          to a slightly perturbed system but the input matrices A, B, D
          and E have not been changed. If SCALE = 0, C and F hold the
          solutions R and L, respectively, to the homogeneous system
          with C = F = 0. Normally, SCALE = 1. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK > = 1.
          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= fla_max(1,2*M*N). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (M+N+6) \n
 * @param[out]	INFO
          INFO is INTEGER \n
            =0: successful exit \n
            <0: If INFO = -i, the i-th argument had an illegal value. \n
            >0: (A, D) and (B, E) have common or close eigenvalues. \n

 *  * */
    template <typename T>
    void tgsyl(char *trans, integer *ijob, integer *m, integer *n, T *a, integer *lda, T *b,
               integer *ldb, T *c, integer *ldc, T *d, integer *ldd, T *e, integer *lde, T *f,
               integer *ldf, T *scale, T *dif, T *work, integer *lwork, integer *iwork,
               integer *info)
    {
        tgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, dif, work,
              lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void tgsyl(char *trans, integer *ijob, integer *m, integer *n, T *a, integer *lda, T *b,
               integer *ldb, T *c, integer *ldc, T *d, integer *ldd, T *e, integer *lde, T *f,
               integer *ldf, Ta *scale, Ta *dif, T *work, integer *lwork, integer *iwork,
               integer *info)
    {
        tgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, dif, work,
              lwork, iwork, info);
    }
    /** @}*/ // end of tgsyl

    /** @defgroup tgsy2 tgsy2
     * @ingroup NYS
     * @{
     */
    /*! @brief TGSY2 solves the generalized Sylvester equation (unblocked algorithm)

 * @details
 * \b Purpose:
    \verbatim
    TGSY2 solves the generalized Sylvester equation:

                A * R - L * B = scale * C                (1)
                D * R - L * E = scale * F,

    using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,
    (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
    N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)
    must be in generalized Schur canonical form, i.e. A, B are upper
    quasi triangular and D, E are upper triangular. The solution (R, L)
    overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor
    chosen to avoid overflow.

    In matrix notation solving equation (1) corresponds to solve
    Z*x = scale*b, where Z is defined as

           Z = [ kron(In, A)  -kron(B**T, Im) ]             (2)
               [ kron(In, D)  -kron(E**T, Im) ],

    Ik is the identity matrix of size k and X**T is the transpose of X.
    kron(X, Y) is the Kronecker product between the matrices X and Y.
    In the process of solving (1), we solve a number of such systems
    where Dim(In), Dim(In) = 1 or 2.

    If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y,
    which is equivalent to solve for R and L in

                A**T * R  + D**T * L   = scale * C           (3)
                R  * B**T + L  * E**T  = scale * -F

    This case is used to compute an estimate of Dif[(A, D), (B, E)] =
    sigma_min(Z) using reverse communication with SLACON.

    STGSY2 also (IJOB >= 1) contributes to the computation in STGSYL
    of an upper bound on the separation between to matrix pairs. Then
    the input (A, D), (B, E) are sub-pencils of the matrix pair in
    STGSYL. See STGSYL for details.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': solve the generalized Sylvester equation (1). \n
          = 'T': solve the 'transposed' system (3). \n
 * @param[in] IJOB
          IJOB is INTEGER \n
          Specifies what kind of functionality to be performed. \n
          = 0: solve (1) only. \n
          = 1: A contribution from this subsystem to a Frobenius
               norm-based estimate of the separation between two matrix
               pairs is computed. (look ahead strategy is used). \n
          = 2: A contribution from this subsystem to a Frobenius
               norm-based estimate of the separation between two matrix
               pairs is computed. (SGECON on sub-systems is used.) \n
          Not referenced if TRANS = 'T'. \n
 * @param[in] M
          M is INTEGER \n
          On entry, M specifies the order of A and D, and the row
          dimension of C, F, R and L. \n
 * @param[in] N
          N is INTEGER \n
          On entry, N specifies the order of B and E, and the column
          dimension of C, F, R and L. \n
 * @param[in] A
          A is REAL array, dimension (LDA, M) \n
          On entry, A contains an upper quasi triangular matrix. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the matrix A. LDA >= fla_max(1, M). \n
 * @param[in] B
          B is REAL array, dimension (LDB, N) \n
          On entry, B contains an upper quasi triangular matrix. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the matrix B. LDB >= fla_max(1, N). \n
 * @param[in,out] C
          C is REAL array, dimension (LDC, N) \n
          On entry, C contains the right-hand-side of the first matrix
          equation in (1). \n
          On exit, if IJOB = 0, C has been overwritten by the
          solution R. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the matrix C. LDC >= fla_max(1, M). \n
 * @param[in] D
          D is REAL array, dimension (LDD, M) \n
          On entry, D contains an upper triangular matrix. \n
 * @param[in] LDD
          LDD is INTEGER \n
          The leading dimension of the matrix D. LDD >= fla_max(1, M). \n
 * @param[in] E
          E is REAL array, dimension (LDE, N) \n
          On entry, E contains an upper triangular matrix. \n
 * @param[in] LDE
          LDE is INTEGER \n
          The leading dimension of the matrix E. LDE >= fla_max(1, N). \n
 * @param[in,out] F
          F is REAL array, dimension (LDF, N) \n
          On entry, F contains the right-hand-side of the second matrix
          equation in (1). \n
          On exit, if IJOB = 0, F has been overwritten by the
          solution L. \n
 * @param[in] LDF
          LDF is INTEGER \n
          The leading dimension of the matrix F. LDF >= fla_max(1, M). \n
 * @param[out] SCALE
          SCALE is REAL \n
          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions
          R and L (C and F on entry) will hold the solutions to a
          slightly perturbed system but the input matrices A, B, D and
          E have not been changed. If SCALE = 0, R and L will hold the
          solutions to the homogeneous system with C = F = 0. Normally,
          SCALE = 1. \n
 * @param[in,out] RDSUM
          RDSUM is REAL \n
          On entry, the sum of squares of computed contributions to
          the Dif-estimate under computation by STGSYL, where the
          scaling factor RDSCAL (see below) has been factored out.
          On exit, the corresponding sum of squares updated with the
          contributions from the current sub-system.
          If TRANS = 'T' RDSUM is not touched. \n
          NOTE: RDSUM only makes sense when STGSY2 is called by STGSYL. \n
 * @param[in,out] RDSCAL
          RDSCAL is REAL \n
          On entry, scaling factor used to prevent overflow in RDSUM. \n
          On exit, RDSCAL is updated w.r.t. the current contributions
          in RDSUM.
          If TRANS = 'T', RDSCAL is not touched. \n
          NOTE: RDSCAL only makes sense when STGSY2 is called by
                STGSYL. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (M+N+2) \n
 * @param[out] PQ
          PQ is INTEGER \n
          On exit, the number of subsystems (of size 2-by-2, 4-by-4 and
          8-by-8) solved by this routine. \n
 * @param[out] INFO
          INFO is INTEGER \n
          On exit, if INFO is set to \n
            =0: Successful exit \n
            <0: If INFO = -i, the i-th argument had an illegal value. \n
            >0: The matrix pairs (A, D) and (B, E) have common or very
                close eigenvalues. \n

 *  * */
    template <typename T>
    void tgsy2(char *trans, integer *ijob, integer *m, integer *n, T *a, integer *lda, T *b,
               integer *ldb, T *c, integer *ldc, T *d, integer *ldd, T *e, integer *lde, T *f,
               integer *ldf, T *scale, T *rdsum, T *rdscal, integer *iwork, integer *pq,
               integer *info)
    {
        tgsy2(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum,
              rdscal, iwork, pq, info);
    }
    template <typename T, typename Ta>
    void tgsy2(char *trans, integer *ijob, integer *m, integer *n, T *a, integer *lda, T *b,
               integer *ldb, T *c, integer *ldc, T *d, integer *ldd, T *e, integer *lde, T *f,
               integer *ldf, Ta *scale, Ta *rdsum, Ta *rdscal, integer *info)
    {
        tgsy2(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum,
              rdscal, info);
    }
    /** @}*/ // end of tgsy2

    /** @defgroup m22 {un, or}m22
     * @ingroup NYS
     * @{
     */
    /*! @brief ORM22 multiplies a general matrix by a banded orthogonal matrix

 * @details
 * \b Purpose:
    \verbatim
    ORM22 overwrites the general real M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'T':      Q**T * C       C * Q**T

    where Q is a real orthogonal matrix of order NQ, with NQ = M if
    SIDE = 'L' and NQ = N if SIDE = 'R'.
    The orthogonal matrix Q processes a 2-by-2 block structure

           [  Q11  Q12  ]
       Q = [            ]
           [  Q21  Q22  ],

    where Q12 is an N1-by-N1 lower triangular matrix and Q21 is an
    N2-by-N2 upper triangular matrix.
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  apply Q (No transpose); \n
          = 'C':  apply Q**T (Conjugate transpose). \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] N1
 * @param[in] N2
          N1 is INTEGER \n
          N2 is INTEGER \n
          The dimension of Q12 and Q21, respectively. N1, N2 >= 0.
          The following requirement must be satisfied: \n
          N1 + N2 = M if SIDE = 'L' and N1 + N2 = N if SIDE = 'R'. \n
 * @param[in] Q
          Q is REAL array, dimension \n
                              (LDQ,M) if SIDE = 'L' \n
                              (LDQ,N) if SIDE = 'R' \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. \n
          LDQ >= fla_max(1,M) if SIDE = 'L'; LDQ >= fla_max(1,N) if SIDE = 'R'. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C. \n
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1)   returns the optimal LWORK. \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If SIDE = 'L', LWORK >= fla_max(1,N); \n
          if SIDE = 'R', LWORK >= fla_max(1,M). \n
          For optimum performance LWORK >= M*N. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void orm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, T *q,
               integer *ldq, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        orm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
    }
    template <typename T>
    void unm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, T *q,
               integer *ldq, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
    }
    /** @}*/ // end of m22

    /** @defgroup lagv2 lagv2
     * @ingroup NYS
     * @{
     */
    /*! @brief LAGV2 computes the Generalized Schur factorization of a real 2-by-2 \n
     matrix pencil (A,B) where B is upper triangular
 * @details
 * \b Purpose:
    \verbatim
    LAGV2 computes the Generalized Schur factorization of a real 2-by-2
    matrix pencil (A,B) where B is upper triangular. This routine
    computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
    SNR such that

    1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0
       types), then

       [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
       [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]

       [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
       [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],

    2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,
       then

       [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
       [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]

       [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
       [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]

       where b11 >= b22 > 0.
    \endverbatim

 * @param[in,out] A \n
          A is REAL array, dimension (LDA, 2)
          On entry, the 2 x 2 matrix A.
          On exit, A is overwritten by the ``A-part'' of the
          generalized Schur form. \n
 * @param[in] LDA
          LDA is INTEGER \n
          THe leading dimension of the array A.  LDA >= 2. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, 2) \n
          On entry, the upper triangular 2 x 2 matrix B.
          On exit, B is overwritten by the ``B-part'' of the
          generalized Schur form. \n
 * @param[in] LDB
          LDB is INTEGER \n
          THe leading dimension of the array B.  LDB >= 2. \n
 * @param[out] ALPHAR
          ALPHAR is REAL array, dimension (2) \n
 * @param[out] ALPHAI
          ALPHAI is REAL array, dimension (2) \n
 * @param[out] BETA
          BETA is REAL array, dimension (2) \n
          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the
          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may
          be zero. \n
 * @param[out] CSL
          CSL is REAL \n
          The cosine of the left rotation matrix. \n
 * @param[out] SNL
          SNL is REAL \n
          The sine of the left rotation matrix. \n
 * @param[out] CSR
          CSR is REAL \n
          The cosine of the right rotation matrix. \n
 * @param[out] SNR
          SNR is REAL \n
          The sine of the right rotation matrix. \n

 * */
    template <typename T>
    void lagv2(T *a, integer *lda, T *b, integer *ldb, T *alphar, T *alphai, T *beta, T *csl,
               T *snl, T *csr, T *snr)
    {
        lagv2(a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr);
    }
    /** @}*/ // end of lagv2

    /** @defgroup tgevc tgevc
     * @ingroup NYS
     * @{
     */
    /*! @brief TGEVC computes some or all of the right and/or left eigenvectors  \n
     of a pair of real matrices
 * @details
 * \b Purpose:
    \verbatim
     TGEVC computes some or all of the right and/or left eigenvectors of
     a pair of real matrices (S,P), where S is a quasi-triangular matrix
     and P is upper triangular.  Matrix pairs of this type are produced by
     the generalized Schur factorization of a matrix pair (A,B):

        A = Q*S*Z**T,  B = Q*P*Z**T

     as computed by SGGHRD + SHGEQZ.

     The right eigenvector x and the left eigenvector y of (S,P)
     corresponding to an eigenvalue w are defined by:

        S*x = w*P*x,  (y**H)*S = w*(y**H)*P,

     where y**H denotes the conjugate tranpose of y.
     The eigenvalues are not input to this routine, but are computed
     directly from the diagonal blocks of S and P.

     This routine returns the matrices X and/or Y of right and left
     eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     where Z and Q are input matrices.
     If Q and Z are the orthogonal factors from the generalized Schur
     factorization of a matrix pair (A,B), then Z*X and Q*Y
     are the matrices of right and left eigenvectors of (A,B).
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'R': compute right eigenvectors only; \n
          = 'L': compute left eigenvectors only; \n
          = 'B': compute both right and left eigenvectors. \n
 * @param[in] HOWMNY
          HOWMNY is CHARACTER*1 \n
          = 'A': compute all right and/or left eigenvectors; \n
          = 'B': compute all right and/or left eigenvectors,
                 backtransformed by the matrices in VR and/or VL; \n
          = 'S': compute selected right and/or left eigenvectors,
                 specified by the logical array SELECT. \n
 * @param[in] SELECT
          SELECT is LOGICAL array, dimension (N) \n
          If HOWMNY='S', SELECT specifies the eigenvectors to be
          computed.  If w(j) is a real eigenvalue, the corresponding
          real eigenvector is computed if SELECT(j) is .TRUE..
          If w(j) and w(j+1) are the real and imaginary parts of a
          complex eigenvalue, the corresponding complex eigenvector
          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,
          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is
          set to .FALSE.. \n
          Not referenced if HOWMNY = 'A' or 'B'. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices S and P.  N >= 0. \n
 * @param[in] S
          S is REAL array, dimension (LDS,N) \n
          The upper quasi-triangular matrix S from a generalized Schur
          factorization, as computed by SHGEQZ. \n
 * @param[in] LDS
          LDS is INTEGER \n
          The leading dimension of array S.  LDS >= fla_max(1,N). \n
 * @param[in] P
          P is REAL array, dimension (LDP,N) \n
          The upper triangular matrix P from a generalized Schur
          factorization, as computed by SHGEQZ.
          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks
          of S must be in positive diagonal form. \n
 * @param[in] LDP
          LDP is INTEGER \n
          The leading dimension of array P.  LDP >= fla_max(1,N). \n
 * @param[in,out] VL
          VL is REAL array, dimension (LDVL,MM) \n
          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
          contain an N-by-N matrix Q (usually the orthogonal matrix Q
          of left Schur vectors returned by SHGEQZ). \n
          On exit, if SIDE = 'L' or 'B', VL contains: \n
          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P); \n
          if HOWMNY = 'B', the matrix Q*Y; \n
          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
                      SELECT, stored consecutively in the columns of
                      VL, in the same order as their eigenvalues. \n
 \n
          A complex eigenvector corresponding to a complex eigenvalue
          is stored in two consecutive columns, the first holding the
          real part, and the second the imaginary part. \n
 \n
          Not referenced if SIDE = 'R'. \n
 * @param[in] LDVL
          LDVL is INTEGER \n
          The leading dimension of array VL.  LDVL >= 1, and if
          SIDE = 'L' or 'B', LDVL >= N. \n
 * @param[in,out] VR
          VR is REAL array, dimension (LDVR,MM) \n
          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
          contain an N-by-N matrix Z (usually the orthogonal matrix Z
          of right Schur vectors returned by SHGEQZ). \n
 \n
          On exit, if SIDE = 'R' or 'B', VR contains: \n
          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P); \n
          if HOWMNY = 'B' or 'b', the matrix Z*X; \n
          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)
                      specified by SELECT, stored consecutively in the
                      columns of VR, in the same order as their
                      eigenvalues. \n
 \n
          A complex eigenvector corresponding to a complex eigenvalue
          is stored in two consecutive columns, the first holding the
          real part and the second the imaginary part. \n
 \n
          Not referenced if SIDE = 'L'. \n
 * @param[in] LDVR
          LDVR is INTEGER \n
          The leading dimension of the array VR.  LDVR >= 1, and if
          SIDE = 'R' or 'B', LDVR >= N. \n
 * @param[in] MM
          MM is INTEGER \n
          The number of columns in the arrays VL and/or VR. MM >= M. \n
 * @param[out] M
          M is INTEGER \n
          The number of columns in the arrays VL and/or VR actually
          used to store the eigenvectors. If HOWMNY = 'A' or 'B', M
          is set to N. Each selected real eigenvector occupies one
          column and each selected complex eigenvector occupies two
          columns. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (6*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex
                eigenvalue. \n

 *  * */
    template <typename T>
    void tgevc(char *side, char *howmny, logical *select, integer *n, T *s, integer *lds, T *p,
               integer *ldp, T *vl, integer *ldvl, T *vr, integer *ldvr, integer *mm, integer *m,
               T *work, integer *info)
    {
        tgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work, info);
    }

    template <typename T, typename Ta>
    void tgevc(char *side, char *howmny, logical *select, integer *n, T *s, integer *lds, T *p,
               integer *ldp, T *vl, integer *ldvl, T *vr, integer *ldvr, integer *mm, integer *m,
               T *work, Ta *rwork, integer *info)
    {
        tgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work, rwork,
              info);
    }
    /** @}*/ // end of tgevc

    /** @defgroup tgexc tgexc
     * @ingroup NYS
     * @{
     */
    /*! @brief TGEXC reorders the generalized real Schur decomposition of a real matrix pair

 * @details
 * \b Purpose:
    \verbatim
     TGEXC reorders the generalized real Schur decomposition of a real
     matrix pair (A,B) using an orthogonal equivalence transformation

                    (A, B) = Q * (A, B) * Z**T,

     so that the diagonal block of (A, B) with row index IFST is moved
     to row ILST.

     (A, B) must be in generalized real Schur canonical form (as returned
     by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     diagonal blocks. B is upper triangular.

     Optionally, the matrices Q and Z of generalized Schur vectors are
     updated.

            Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
            Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T
    \endverbatim

 * @param[in] WANTQ
          WANTQ is LOGICAL \n
          .TRUE. : update the left transformation matrix Q; \n
          .FALSE.: do not update Q. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          .TRUE. : update the right transformation matrix Z; \n
          .FALSE.: do not update Z. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the matrix A in generalized real Schur canonical
          form. \n
          On exit, the updated matrix A, again in generalized
          real Schur canonical form. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the matrix B in generalized real Schur canonical
          form (A,B). \n
          On exit, the updated matrix B, again in generalized
          real Schur canonical form (A,B). \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,N). \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.
          On exit, the updated matrix Q. \n
          If WANTQ = .FALSE., Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. LDQ >= 1.
          If WANTQ = .TRUE., LDQ >= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,N) \n
          On entry, if WANTZ = .TRUE., the orthogonal matrix Z.
          On exit, the updated matrix Z. \n
          If WANTZ = .FALSE., Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z. LDZ >= 1.
          If WANTZ = .TRUE., LDZ >= N. \n
 * @param[in,out] IFST
          IFST is INTEGER \n
 * @param[in,out] ILST
          ILST is INTEGER \n
          Specify the reordering of the diagonal blocks of (A, B).
          The block with row index IFST is moved to row ILST, by a
          sequence of swapping between adjacent blocks.
          On exit, if IFST pointed on entry to the second row of
          a 2-by-2 block, it is changed to point to the first row;
          ILST always points to the first row of the block in its
          final position (which may differ from its input value by
          +1 or -1). 1 <= IFST, ILST <= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.
          LWORK >= 1 when N <= 1, otherwise LWORK >= 4*N + 16. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
           =0:  successful exit. \n
           <0:  if INFO = -i, the i-th argument had an illegal value. \n
           =1:  The transformed matrix pair (A, B) would be too far
                from generalized Schur form; the problem is ill-
                conditioned. (A, B) may have been partially reordered,
                and ILST points to the first row of the current
                position of the block being moved. \n

 *  * */
    template <typename T>
    void tgexc(logical *wantq, logical *wantz, integer *n, T *a, integer *lda, T *b, integer *ldb,
               T *q, integer *ldq, T *z, integer *ldz, integer *ifst, integer *ilst, integer *info)
    {
        tgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info);
    }

    template <typename T, typename Ta>
    void tgexc(logical *wantq, logical *wantz, integer *n, T *a, integer *lda, T *b, integer *ldb,
               T *q, integer *ldq, T *z, integer *ldz, integer *ifst, integer *ilst, Ta *work,
               integer *lwork, integer *info)
    {
        tgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info);
    }
    /** @}*/ // end of tgexc

    /** @defgroup tgex2 tgex2
     * @ingroup NYS
     * @{
     */
    /*! @brief TGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an
 orthogonal equivalence transformation

 * @details
 * \b Purpose:
    \verbatim
    TGEX2 swaps adjacent diagonal blocks (A11, B11) and (A22, B22)
    of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
    (A, B) by an orthogonal equivalence transformation.

    (A, B) must be in generalized real Schur canonical form (as   returned
    by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
    diagonal blocks. B is upper triangular.

    Optionally, the matrices Q and Z of generalized Schur vectors are
    updated.

           Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
           Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T
    \endverbatim

 * @param[in] WANTQ
          WANTQ is LOGICAL \n
          .TRUE. : update the left transformation matrix Q;
          .FALSE.: do not update Q. \n
 * @param[in] WANTZ
          WANTZ is LOGICAL \n
          .TRUE. : update the right transformation matrix Z;
          .FALSE.: do not update Z. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B. N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the matrix A in the pair (A, B).
          On exit, the updated matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the matrix B in the pair (A, B).
          On exit, the updated matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,N). \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.
          On exit, the updated matrix Q.
          Not referenced if WANTQ = .FALSE.. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. LDQ >= 1.
          If WANTQ = .TRUE., LDQ >= N. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ,N) \n
          On entry, if WANTZ =.TRUE., the orthogonal matrix Z.
          On exit, the updated matrix Z.
          Not referenced if WANTZ = .FALSE.. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z. LDZ >= 1.
          If WANTZ = .TRUE., LDZ >= N. \n
 * @param[in] J1
          J1 is INTEGER \n
          The index to the first block (A11, B11). 1 <= J1 <= N. \n
 * @param[in] N1
          N1 is INTEGER \n
          The order of the first block (A11, B11). N1 = 0, 1 or 2. \n
 * @param[in] N2
          N2 is INTEGER \n
          The order of the second block (A22, B22). N2 = 0, 1 or 2. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)). \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.
          LWORK >=  MAX(N*(N2+N1), (N2+N1)*(N2+N1)*2) \n
 * @param[out] INFO
          INFO is INTEGER \n
            =0: Successful exit \n
            >0: If INFO = 1, the transformed matrix (A, B) would be
                too far from generalized Schur form; the blocks are
                not swapped and (A, B) and (Q, Z) are unchanged.
                The problem of swapping is too ill-conditioned. \n
            <0: If INFO = -16: LWORK is too small. Appropriate value
                for LWORK is   returned in WORK(1).  \n

 *  * */
    template <typename T>
    void tgex2(logical *wantq, logical *wantz, integer *n, T *a, integer *lda, T *b, integer *ldb,
               T *q, integer *ldq, T *z, integer *ldz, integer *j1, integer *n1, integer *n2,
               T *work, integer *lwork, integer *info)
    {
        tgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info);
    }
    /** @}*/ // end of tgex2
    /**@}*/ // end of NYS

    /** @defgroup SYM Symmetric eigenvalues
     * @ingroup Eig
     * @{
     */
    /** @defgroup syev syev
     * @ingroup SYM
     * @{
     */
    /*! @brief Eigenvalue decomposition (QR algorithm)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Eigenvalue decomposition (QR algorithm).
        Computation of all eigenvalues and, optionally, eigenvectors of a real symmetric matrix a.
    \endverbatim

    * @param[in] jobz
              jobz is char* \n
              = 'N':  Compute eigenvalues only; \n
              = 'V':  Compute eigenvalues and eigenvectors. \n
    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda, n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the
              leading n-by-n upper triangular part of a contains the
              upper triangular part of the matrix a.  If uplo = 'L',
              the leading n-by-n lower triangular part of a contains
              the lower triangular part of the matrix a. \n
              On exit, if jobz = 'V', then if info = 0, A contains the
              orthonormal eigenvectors of the matrix a. \n
              If jobz = 'N', then on exit the lower triangle (if uplo='L')
              or the upper triangle (if uplo='U') of A, including the
              diagonal, is destroyed. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] w
              w is float/double array, dimension (n) \n
              If info = 0, the eigenvalues in ascending order. \n
    * @param[out]	WORK
              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The length of the array WORK.  LWORK >= fla_max(1,3*N-1). \n
              For optimal efficiency, LWORK >= (NB+2)*N,
              where NB is the blocksize for DSYTRD returned by ILAENV. \n

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  if INFO = i, the algorithm failed to converge; i
                    off-diagonal elements of an intermediate tridiagonal
                    form did not converge to zero. \n

    *     *  */
    template <typename T>
    void syev(char *jobz, char *uplo, integer *n, T *a, integer *lda, T *w, T *work, integer *lwork,
              integer *info)
    {
        syev(jobz, uplo, n, a, lda, w, work, lwork, info);
    }
    /** @}*/ // end of syev

    /** @defgroup syevd syevd
     * @ingroup SYM
     * @{
     */
    /*! @brief Eigenvalue decomposition (divide-and-conquer)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of the eigenvalues and, optionally, the left and/or right eigenvectors for SY
        matrices. SSYEVD computes all eigenvalues and, optionally, eigenvectors of a real
        symmetric matrix a. If eigenvectors are desired, it uses a divide and conquer algorithm.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP, Cray
        C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.

        Because of large use of BLAS of level 3, SSYEVD needs N**2 more workspace than SSYEVX.
    \endverbatim

    * @param[in] jobz
              jobz is char* \n
              = 'N':  Compute eigenvalues only; \n
              = 'V':  Compute eigenvalues and eigenvectors. \n
    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda, n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the
              leading n-by-n upper triangular part of a contains the
              upper triangular part of the matrix a.  If uplo = 'L',
              the leading n-by-n lower triangular part of a contains
              the lower triangular part of the matrix a. \n
              On exit, if jobz = 'V', then if info = 0, A contains the
              orthonormal eigenvectors of the matrix a. \n
              If jobz = 'N', then on exit the lower triangle (if uplo='L')
              or the upper triangle (if uplo='U') of A, including the
              diagonal, is destroyed. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] w
              w is float/double array, dimension (n) \n
              If info = 0, the eigenvalues in ascending order. \n
    * @param[out]	WORK
              WORK is DOUBLE PRECISION array,
                                             dimension (LWORK) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              If N <= 1,               LWORK must be at least 1. \n
              If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1. \n
              If JOBZ = 'V' and N > 1, LWORK must be at least
                                                    1 + 6*N + 2*N**2. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal sizes of the WORK and IWORK
              arrays, returns these values as the first entries of the WORK
              and IWORK arrays, and no error message related to LWORK or
              LIWORK is issued by XERBLA. \n
    * @param[out]	IWORK
              IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
              On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
    * @param[in]	LIWORK
              LIWORK is INTEGER \n
              The dimension of the array IWORK. \n
              If N <= 1,                LIWORK must be at least 1. \n
              If JOBZ  = 'N' and N > 1, LIWORK must be at least 1. \n
              If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. \n
 \n
              If LIWORK = -1, then a workspace query is assumed; the
              routine only calculates the optimal sizes of the WORK and
              IWORK arrays, returns these values as the first entries of
              the WORK and IWORK arrays, and no error message related to
              LWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
                    to converge; i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero;
                    if INFO = i and JOBZ = 'V', then the algorithm failed
                    to compute an eigenvalue while working on the submatrix
                    lying in rows and columns INFO/(N+1) through
                    mod(INFO,N+1). \n

    *     *  */
    template <typename T>
    void syevd(char *jobz, char *uplo, integer *n, T *a, integer *lda, T *w, T *work,
               integer *lwork, integer *iwork, integer *liwork, integer *info)
    {
        syevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
    }
    /** @}*/ // end of syevd

    /** @defgroup syevr syevr
     * @ingroup SYM
     * @{
     */
    /*! @brief Hermitian eigenvalue decomposition (MRRR)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Hermitian eigenvalue decomposition (MRRR).
        Computation of eigenvalues and, optionally, the left and/or right eigenvectors for SY
        matrices

        SSYEVR computes selected eigenvalues and, optionally, eigenvectors of a real symmetric
        matrix a. Eigenvalues and eigenvectors can be selected by specifying either a range of
        values or a range of indices for the desired eigenvalues.

        SSYEVR first reduces the matrix a to tridiagonal form T with a call to SSYTRD.  Then,
        whenever possible, SSYEVR calls SSTEMR to compute the eigenspectrum using Relatively
        Robust Representations. SSTEMR computes eigenvalues by the dqds algorithm, while orthogonal
        eigenvectors are computed from various "good" L D L^T representations (also known as
        Relatively Robust Representations). Gram-Schmidt orthogonalization is avoided as far as
        possible. More specifically, the various steps of the algorithm are as follows.

        For each unreduced block (submatrix) of T,
        (a) Compute T - sigma I  = L D L^T, so that L and D
        define all the wanted eigenvalues to high relative accuracy.
        This means that small relative changes in the entries of D and L
        cause only small relative changes in the eigenvalues and
        eigenvectors. The standard (unfactored) representation of the
        tridiagonal matrix T does not have this property in general.
        (b) Compute the eigenvalues to suitable accuracy.
        If the eigenvectors are desired, the algorithm attains full
        accuracy of the computed eigenvalues only right before
        the corresponding vectors have to be computed, see steps c) and d).
        (c) For each cluster of close eigenvalues, select a new
        shift close to the cluster, find a new factorization, and refine
        the shifted eigenvalues to suitable accuracy.
        (d) For each eigenvalue with a large enough relative separation compute
        the corresponding eigenvector by forming a rank revealing twisted
        factorization. Go back to (c) for any clusters that remain.

        The desired accuracy of the output can be specified by the input parameter abstol.

        For more details, see SSTEMR's documentation and:
        - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations to compute
        orthogonal eigenvectors of symmetric tridiagonal matrices," Linear Algebra and its
        Applications, 387(1), pp. 1-28, August 2004.
        - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and   Relative Gaps,"
        SIAM Journal on Matrix Analysis and Applications, Vol. 25,   2004.
        Also LAPACK Working Note 154.
        - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric tridiagonal
        eigenvalue/eigenvector problem",   Computer Science Division Technical Report No.
        UCB/CSD-97-971, UC Berkeley, May 1997.


        Note 1 : SSYEVR calls SSTEMR when the full spectrum is requested on machines which conform
        to the ieee-754 floating point standard. SSYEVR calls SSTEBZ and SSTEIN on non-ieee
        machines and when partial spectrum requests are made.

        Normal execution of SSTEMR may create NaNs and infinities and hence may abort due to a
        floating point exception in environments which do not handle NaNs and infinities in the
        ieee standard default manner.
    \endverbatim

    * @param[in] jobz
              jobz is char* \n
              = 'N':  Compute eigenvalues only; \n
              = 'V':  Compute eigenvalues and eigenvectors. \n
    * @param[in] range
              range is char* \n
              = 'A': all eigenvalues will be found. \n
              = 'V': all eigenvalues in the half-open interval (vl,vu]
              will be found. \n
              = 'I': the il-th through iu-th eigenvalues will be found. \n
              For range = 'V' or 'I' and iu - il < N - 1, SSTEBZ and
              SSTEIN are called \n
    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda, n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the
              leading n-by-n upper triangular part of a contains the
              upper triangular part of the matrix a.  If uplo = 'L',
              the leading n-by-n lower triangular part of a contains
              the lower triangular part of the matrix a. \n
              On exit, the lower triangle (if uplo='L') or the upper
              triangle (if uplo='U') of A, including the diagonal, is
              destroyed. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[in] vl
              vl is float/double* \n
              If range='V', the lower bound of the interval to
              be searched for eigenvalues. vl < vu. \n
              Not referenced if range = 'A' or 'I'. \n
    * @param[in] vu
              vu is float/double* \n
              If range='V', the upper bound of the interval to
              be searched for eigenvalues. vl < vu. \n
              Not referenced if range = 'A' or 'I'. \n
    * @param[in] il
              il is integer* \n
              If range='I', the index of the
              smallest eigenvalue to be returned. \n
              1 <= il <= iu <= N, if N > 0; il = 1 and iu = 0 if N = 0. \n
              Not referenced if range = 'A' or 'V'. \n
    * @param[in] iu
              iu is integer* \n
              If range='I', the index of the
              largest eigenvalue to be returned. \n
              1 <= il <= iu <= N, if N > 0; il = 1 and iu = 0 if N = 0. \n
              Not referenced if range = 'A' or 'V'. \n
    * @param[in] abstol
              abstol is float/double* \n
              The absolute error tolerance for the eigenvalues. \n
              An approximate eigenvalue is accepted as converged
              when it is determined to lie in an interval [a,b]
              of width less than or equal to \n \n
              abstol + EPS *   fla_max( |a|,|b|) , \n \n
              where EPS is the machine precision.  If abstol is less than
              or equal to zero, then  EPS*|T|  will be used in its place,
              where |T| is the 1-norm of the tridiagonal matrix obtained
              by reducing A to tridiagonal form. \n \n
              See "Computing Small Singular Values of Bidiagonal Matrices
              with Guaranteed High Relative Accuracy," by Demmel and
              Kahan, LAPACK Working Note #3. \n \n
              If high relative accuracy is important, set abstol to
              SLAMCH( 'Safe minimum').  Doing so will guarantee that
              eigenvalues are computed to high relative accuracy when
              possible in future releases.  The current code does not
              make any guarantees about high relative accuracy, but
              future releases will. See J. Barlow and J. Demmel,
              "Computing Accurate Eigensystems of Scaled Diagonally
              Dominant Matrices", LAPACK Working Note #7, for a discussion
              of which matrices define their eigenvalues to high relative
              accuracy. \n
    * @param[out] m
              m is integer* \n
              The total number of eigenvalues found.  0 <= m <= n. \n
              If range = 'A', m = n, and if range = 'I', m = iu-il+1. \n
    * @param[out] w
              w is float/double/COMPLEX/COMPLEX*16 array, dimension (n) \n
              The first m elements contain the selected eigenvalues in
              ascending order. \n
    * @param[out] z
              z is float/double array, dimension (ldz, fla_max(1,m)) \n
              If jobz = 'V', then if info = 0, the first M columns of Z
              contain the orthonormal eigenvectors of the matrix a
              corresponding to the selected eigenvalues, with the i-th
              column of Z holding the eigenvector associated with W(i). \n
              If jobz = 'N', then Z is not referenced. \n
              Note: the user must ensure that at least fla_max(1,m) columns are
              supplied in the array Z; if range = 'V', the exact value of m
              is not known in advance and an upper bound must be used.
              Supplying n columns is always safe. \n
    * @param[in] ldz
              ldz is integer* \n
              The leading dimension of the array Z.  ldz >= 1, and if
              jobz = 'V', ldz >= fla_max(1,n). \n
    * @param[out] isuppz
              isuppz is integer array, dimension ( 2*max(1,m)) \n
              The support of the eigenvectors in Z, i.e., the indices
              indicating the nonzero elements in Z. The i-th eigenvector
              is nonzero only in elements isuppz( 2*i-1) through
              isuppz( 2*i). This is an output of SSTEMR (tridiagonal
              matrix). The support of the eigenvectors of A is typically
              1:N because of the orthogonal transformations applied by SORMTR. \n
              Implemented only for range = 'A' or 'I' and iu - il = N - 1 \n
    * @param[out]	WORK
              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK.  LWORK >= fla_max(1,26*N).
              For optimal efficiency, LWORK >= (NB+6)*N,
              where NB is the max of the blocksize for DSYTRD and DORMTR
              returned by ILAENV. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	IWORK
              IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
              On exit, if INFO = 0, IWORK(1) returns the optimal LWORK. \n
    * @param[in]	LIWORK
              LIWORK is INTEGER \n
              The dimension of the array IWORK.  LIWORK >= fla_max(1,10*N). \n
 \n
              If LIWORK = -1, then a workspace query is assumed; the
              routine only calculates the optimal size of the IWORK array,
              returns this value as the first entry of the IWORK array, and
              no error message related to LIWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  Internal error \n

    *     *  */
    template <typename T>
    void syevr(char *jobz, char *range, char *uplo, integer *n, T *a, integer *lda, T *vl, T *vu,
               integer *il, integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz,
               integer *isuppz, T *work, integer *lwork, integer *iwork, integer *liwork,
               integer *info)
    {
        syevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
              lwork, iwork, liwork, info);
    }
    /** @}*/ // end of syevr
    /** @defgroup syevx syevx
     * @ingroup SYM
     * @{
     */
    /*! @brief SYEVX  computes the eigenvalues and, optionally, the left \n
     and/or right eigenvectors for SY matrices
 * @details
 * \b Purpose:
    \verbatim
     SYEVX computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
     selected by specifying either a range of values or a range of indices
     for the desired eigenvalues.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A. \n
          On exit, the lower triangle (if UPLO='L') or the upper
          triangle (if UPLO='U') of A, including the diagonal, is
          destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*SLAMCH('S'). \n
 \n
          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          On normal exit, the first M elements contain the selected
          eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          If an eigenvector fails to converge, then that column of Z
          contains the latest approximation to the eigenvector, and the
          index of the eigenvector is returned in IFAIL. \n
          If JOBZ = 'N', then Z is not referenced. \n
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of the array WORK.  LWORK >= 1, when N <= 1;
          otherwise 8*N. \n
          For optimal efficiency, LWORK >= (NB+3)*N,
          where NB is the max of the blocksize for SSYTRD and SORMTR
          returned by ILAENV. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (5*N) \n
 * @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (N) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, then i eigenvectors failed to converge.
                Their indices are stored in array IFAIL. \n

 *  * */
    template <typename T>
    void syevx(char *jobz, char *range, char *uplo, integer *n, T *a, integer *lda, T *vl, T *vu,
               integer *il, integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz, T *work,
               integer *lwork, integer *iwork, integer *ifail, integer *info)
    {
        syevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
              iwork, ifail, info);
    }
    /** @}*/ // end of syevx

    /** @defgroup syevx_2stage syevx_2stage
     * @ingroup SYM
     * @{
     */
    /*! @brief SYEVX_2STAGE  computes the eigenvalues and, optionally, the left \n
     and/or right eigenvectors for SY matrices
 * @details
 * \b Purpose:
    \verbatim
     SYEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric matrix A using the 2stage technique for
     the reduction to tridiagonal.  Eigenvalues and eigenvectors can be
     selected by specifying either a range of values or a range of indices
     for the desired eigenvalues.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors.
                  Not available in this release. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A. \n
          On exit, the lower triangle (if UPLO='L') or the upper
          triangle (if UPLO='U') of A, including the diagonal, is
          destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*SLAMCH('S'). \n
 \n
          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          On normal exit, the first M elements contain the selected
          eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          If an eigenvector fails to converge, then that column of Z
          contains the latest approximation to the eigenvector, and the
          index of the eigenvector is returned in IFAIL. \n
          If JOBZ = 'N', then Z is not referenced.
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of the array WORK. LWORK >= 1, when N <= 1; \n
          otherwise   \n
          If JOBZ = 'N' and N > 1, LWORK must be queried.
                                   LWORK = MAX(1, 8*N, dimension) where
                                   dimension = fla_max(stage1,stage2) + (KD+1)*N + 3*N
                                             = N*KD + N*max(KD+1,FACTOPTNB)
                                               + fla_max(2*KD*KD, KD*NTHREADS)
                                               + (KD+1)*N + 3*N
                                   where KD is the blocking size of the reduction,
                                   FACTOPTNB is the blocking used by the QR or LQ
                                   algorithm, usually FACTOPTNB=128 is a good choice
                                   NTHREADS is the number of threads used when
                                   openMP compilation is enabled, otherwise =1. \n
          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (5*N) \n
 * @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (N) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, then i eigenvectors failed to converge.
                Their indices are stored in array IFAIL. \n

 *  * */
    template <typename T>
    void syevx_2stage(char *jobz, char *range, char *uplo, integer *n, T *a, integer *lda, T *vl,
                      T *vu, integer *il, integer *iu, T *abstol, integer *m, T *w, T *z,
                      integer *ldz, T *work, integer *lwork, integer *iwork, integer *ifail,
                      integer *info)
    {
        syevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work,
                     lwork, iwork, ifail, info);
    }
    /** @}*/ // end of syevx_@stage

    /** @defgroup syev_2stage syev_2stage
     * @ingroup SYM
     * @{
     */
    /*! @brief SYEV_2STAGE computes the eigenvalues and, optionally, the left \n
     and/or right eigenvectors for SY matrices
 * @details
 * \b Purpose:
    \verbatim
     SYEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of a
     real symmetric matrix A using the 2stage technique for
     the reduction to tridiagonal.
    \endverbatim

  * @param[in] JOBZ
           JOBZ is CHARACTER*1 \n
           = 'N':  Compute eigenvalues only; \n
           = 'V':  Compute eigenvalues and eigenvectors.
                   Not available in this release. \n
  * @param[in] UPLO
           UPLO is CHARACTER*1 \n
           = 'U':  Upper triangle of A is stored; \n
           = 'L':  Lower triangle of A is stored. \n
  * @param[in] N
           N is INTEGER \n
           The order of the matrix A.  N >= 0. \n
  * @param[in,out] A
           A is REAL array, dimension (LDA, N) \n
           On entry, the symmetric matrix A.  If UPLO = 'U', the
           leading N-by-N upper triangular part of A contains the
           upper triangular part of the matrix A.  If UPLO = 'L',
           the leading N-by-N lower triangular part of A contains
           the lower triangular part of the matrix A. \n
           On exit, if JOBZ = 'V', then if INFO = 0, A contains the
           orthonormal eigenvectors of the matrix A.
           If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
           or the upper triangle (if UPLO='U') of A, including the
           diagonal, is destroyed. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,N). \n
  * @param[out] W
           W is REAL array, dimension (N) \n
           If INFO = 0, the eigenvalues in ascending order. \n
  * @param[out]	WORK
          WORK is REAL array, dimension LWORK \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
  * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of the array WORK. LWORK >= 1, when N <= 1; \n
          otherwise   \n
          If JOBZ = 'N' and N > 1, LWORK must be queried.
                                   LWORK = MAX(1, dimension) where \n
                                   dimension = fla_max(stage1,stage2) + (KD+1)*N + 2*N \n
                                             = N*KD + N*max(KD+1,FACTOPTNB)
                                               + fla_max(2*KD*KD, KD*NTHREADS)
                                               + (KD+1)*N + 2*N \n
                                   where KD is the blocking size of the reduction,
                                   FACTOPTNB is the blocking used by the QR or LQ
                                   algorithm, usually FACTOPTNB=128 is a good choice
                                   NTHREADS is the number of threads used when
                                   openMP compilation is enabled, otherwise =1. \n
          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
  * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the algorithm failed to converge; i
                off-diagonal elements of an intermediate tridiagonal
                form did not converge to zero. \n

 *  * */
    template <typename T>
    void syev_2stage(char *jobz, char *uplo, integer *n, T *a, integer *lda, T *w, T *work,
                     integer *lwork, integer *info)
    {
        syev_2stage(jobz, uplo, n, a, lda, w, work, lwork, info);
    }
    /** @}*/ // end of syev2stage

    /** @defgroup syevr_2stage syevr_2stage
     * @ingroup SYM
     * @{
     */
    /*! @brief SYEVR_2STAGE computes the eigenvalues and, optionally, the left \n
     and/or right eigenvectors for SY matrices
 * @details
 * \b Purpose:
    \verbatim
     SYEVR_2STAGE computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric matrix A using the 2stage technique for
     the reduction to tridiagonal.  Eigenvalues and eigenvectors can be
     selected by specifying either a range of values or a range of
     indices for the desired eigenvalues.

     SYEVR_2STAGE first reduces the matrix A to tridiagonal form T with a call
     to SYTRD.  Then, whenever possible, SYEVR_2STAGE calls STEMR to compute
     the eigenspectrum using Relatively Robust Representations.  STEMR
     computes eigenvalues by the dqds algorithm, while orthogonal
     eigenvectors are computed from various "good" L D L^T representations
     (also known as Relatively Robust Representations). Gram-Schmidt
     orthogonalization is avoided as far as possible. More specifically,
     the various steps of the algorithm are as follows.

     For each unreduced block (submatrix) of T,
        (a) Compute T - sigma I  = L D L^T, so that L and D
            define all the wanted eigenvalues to high relative accuracy.
            This means that small relative changes in the entries of D and L
            cause only small relative changes in the eigenvalues and
            eigenvectors. The standard (unfactored) representation of the
            tridiagonal matrix T does not have this property in general.
        (b) Compute the eigenvalues to suitable accuracy.
            If the eigenvectors are desired, the algorithm attains full
            accuracy of the computed eigenvalues only right before
            the corresponding vectors have to be computed, see steps c) and d).
        (c) For each cluster of close eigenvalues, select a new
            shift close to the cluster, find a new factorization, and refine
            the shifted eigenvalues to suitable accuracy.
        (d) For each eigenvalue with a large enough relative separation compute
            the corresponding eigenvector by forming a rank revealing twisted
            factorization. Go back to (c) for any clusters that remain.

     The desired accuracy of the output can be specified by the input
     parameter ABSTOL.

     For more details, see SSTEMR's documentation and:
     - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
       to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
       Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
     - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
       Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
       2004.  Also LAPACK Working Note 154.
     - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
       tridiagonal eigenvalue/eigenvector problem",
       Computer Science Division Technical Report No. UCB/CSD-97-971,
       UC Berkeley, May 1997.


     Note 1 : SSYEVR_2STAGE calls SSTEMR when the full spectrum is requested
     on machines which conform to the ieee-754 floating point standard.
     SSYEVR_2STAGE calls SSTEBZ and SSTEIN on non-ieee machines and
     when partial spectrum requests are made.

     Normal execution of SSTEMR may create NaNs and infinities and
     hence may abort due to a floating point exception in environments
     which do not handle NaNs and infinities in the ieee standard default
     manner.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors.
                  Not available in this release. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found.
          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and
          SSTEIN are called \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A. \n
          On exit, the lower triangle (if UPLO='L') or the upper
          triangle (if UPLO='U') of A, including the diagonal, is
          destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form. \n
 \n
          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3. \n
 \n
          If high relative accuracy is important, set ABSTOL to
          SLAMCH( 'Safe minimum').  Doing so will guarantee that
          eigenvalues are computed to high relative accuracy when
          possible in future releases.  The current code does not
          make any guarantees about high relative accuracy, but
          future releases will. See J. Barlow and J. Demmel,
          "Computing Accurate Eigensystems of Scaled Diagonally
          Dominant Matrices", LAPACK Working Note #7, for a discussion
          of which matrices define their eigenvalues to high relative
          accuracy. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first M elements contain the selected eigenvalues in
          ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i). \n
          If JOBZ = 'N', then Z is not referenced. \n
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used.
          Supplying N columns is always safe. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out] ISUPPZ
          ISUPPZ is INTEGER array, dimension ( 2*max(1,M)) \n
          The support of the eigenvectors in Z, i.e., the indices
          indicating the nonzero elements in Z. The i-th eigenvector
          is nonzero only in elements ISUPPZ( 2*i-1) through
          ISUPPZ( 2*i). This is an output of SSTEMR (tridiagonal
          matrix). The support of the eigenvectors of A is typically
          1:N because of the orthogonal transformations applied by SORMTR.
          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If JOBZ = 'N' and N > 1, LWORK must be queried.
                                   LWORK = MAX(1, 26*N, dimension) where
                                   dimension = fla_max(stage1,stage2) + (KD+1)*N + 5*N
                                             = N*KD + N*max(KD+1,FACTOPTNB)
                                               + fla_max(2*KD*KD, KD*NTHREADS)
                                               + (KD+1)*N + 5*N
                                   where KD is the blocking size of the reduction,
                                   FACTOPTNB is the blocking used by the QR or LQ
                                   algorithm, usually FACTOPTNB=128 is a good choice
                                   NTHREADS is the number of threads used when
                                   openMP compilation is enabled, otherwise =1. \n
          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK.  LIWORK >= fla_max(1,10*N). \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the IWORK array,
          returns this value as the first entry of the IWORK array, and
          no error message related to LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  Internal error \n

 *  * */
    template <typename T>
    void syevr_2stage(char *jobz, char *range, char *uplo, integer *n, T *a, integer *lda, T *vl,
                      T *vu, integer *il, integer *iu, T *abstol, integer *m, T *w, T *z,
                      integer *ldz, integer *isuppz, T *work, integer *lwork, integer *iwork,
                      integer *liwork, integer *info)
    {
        syevr_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz,
                     work, lwork, iwork, liwork, info);
    }
    /** @}*/ // end of

    /** @defgroup spev spev
     * @ingroup SYM
     * @{
     */
    /*! @brief SSPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     SSPEV computes all the eigenvalues and, optionally, eigenvectors of a
     real symmetric matrix A in packed storage.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
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
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, AP is overwritten by values generated during the
          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
          and first superdiagonal of the tridiagonal matrix T overwrite
          the corresponding elements of A, and if UPLO = 'L', the
          diagonal and first subdiagonal of T overwrite the
          corresponding elements of A. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
          eigenvectors of the matrix A, with the i-th column of Z
          holding the eigenvector associated with W(i). \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, the algorithm failed to converge; i
                off-diagonal elements of an intermediate tridiagonal
                form did not converge to zero. \n

 *  * */
    template <typename T>
    void spev(char *jobz, char *uplo, integer *n, T *ap, T *w, T *z, integer *ldq, T *work,
              integer *info)
    {
        spev(jobz, uplo, n, ap, w, z, ldq, work, info);
    }
    template <typename T, typename Ta>
    void hpev(char *jobz, char *uplo, integer *n, T *ap, Ta *w, T *z, integer *ldq, T *work,
              Ta *rwork, integer *info)
    {
        hpev(jobz, uplo, n, ap, w, z, ldq, work, rwork, info);
    }
    /** @}*/ // end of spev

    /** @defgroup spevd spevd
     * @ingroup SYM
     * @{
     */
    /*! @brief SPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors
 for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     SPEVD computes all the eigenvalues and, optionally, eigenvectors
     of a real symmetric matrix A in packed storage. If eigenvectors are
     desired, it uses a divide and conquer algorithm.

     The divide and conquer algorithm makes very mild assumptions about
     floating point arithmetic. It will work on machines with a guard
     digit in add/subtract, or on those binary machines without guard
     digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     Cray-2. It could conceivably fail on hexadecimal or decimal machines
     without guard digits, but we know of none.
    \endverbatim

 * @param[in] JOBZ \n
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
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
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, AP is overwritten by values generated during the
          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
          and first superdiagonal of the tridiagonal matrix T overwrite
          the corresponding elements of A, and if UPLO = 'L', the
          diagonal and first subdiagonal of T overwrite the
          corresponding elements of A. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
          eigenvectors of the matrix A, with the i-th column of Z
          holding the eigenvector associated with W(i). \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the required LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If N <= 1,               LWORK must be at least 1. \n
          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N. \n
          If JOBZ = 'V' and N > 1, LWORK must be at least
                                                 1 + 6*N + N**2. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the required sizes of the WORK and IWORK
          arrays, returns these values as the first entries of the WORK
          and IWORK arrays, and no error message related to LWORK or
          LIWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the required LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1. \n
          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the required sizes of the WORK and
          IWORK arrays, returns these values as the first entries of
          the WORK and IWORK arrays, and no error message related to
          LWORK or LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, the algorithm failed to converge; i
                off-diagonal elements of an intermediate tridiagonal
                form did not converge to zero. \n

 *  * */
    template <typename T>
    void spevd(char *jobz, char *uplo, integer *n, T *ap, T *w, T *z, integer *ldz, T *work,
               integer *lwork, integer *iwork, integer *liwork, integer *info)
    {
        spevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
    }
    template <typename T, typename Ta>
    void hpevd(char *jobz, char *uplo, integer *n, T *ap, Ta *w, T *z, integer *ldz, T *work,
               integer *lwork, Ta *rwork, integer *lrwork, integer *iwork, integer *liwork,
               integer *info)
    {
        hpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
    }

    /** @}*/ // end of spevd

    /** @defgroup spevx spevx
     * @ingroup SYM
     * @{
     */
    /*! @brief SPEVX computes the eigenvalues and, optionally, the left and/or right
 eigenvectors for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     SPEVX computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
     can be selected by specifying either a range of values or a range of
     indices for the desired eigenvalues.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found; \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found; \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
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
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, AP is overwritten by values generated during the
          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
          and first superdiagonal of the tridiagonal matrix T overwrite
          the corresponding elements of A, and if UPLO = 'L', the
          diagonal and first subdiagonal of T overwrite the
          corresponding elements of A. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing AP to tridiagonal form. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*SLAMCH('S'). \n
 \n
          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the selected eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          If an eigenvector fails to converge, then that column of Z
          contains the latest approximation to the eigenvector, and the
          index of the eigenvector is returned in IFAIL. \n
          If JOBZ = 'N', then Z is not referenced. \n
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (8*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (5*N) \n
 * @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (N) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, then i eigenvectors failed to converge.
                Their indices are stored in array IFAIL. \n

 *  * */
    template <typename T>
    void spevx(char *jobz, char *range, char *uplo, integer *n, T *ap, T *vl, T *vu, integer *il,
               integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz, T *work,
               integer *iwork, integer *ifail, integer *info)
    {
        spevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail,
              info);
    }
    template <typename T, typename Ta>
    void spevx(char *jobz, char *range, char *uplo, integer *n, T *ap, Ta *vl, Ta *vu, integer *il,
               integer *iu, Ta *abstol, integer *m, Ta *w, T *z, integer *ldz, T *work, Ta *rwork,
               integer *iwork, integer *ifail, integer *info)
    {
        spevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork,
              ifail, info);
    }
    /** @}*/ // end of spevx

    /** @defgroup stev stev
     * @ingroup SYM
     * @{
     */
    /*! @brief STEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     STEV computes all eigenvalues and, optionally, eigenvectors of a
     real symmetric tridiagonal matrix A.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the n diagonal elements of the tridiagonal matrix
          A.
          On exit, if INFO = 0, the eigenvalues in ascending order. \n
 * @param[in,out] E
          E is REAL array, dimension (N-1) \n
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix A, stored in elements 1 to N-1 of E.
          On exit, the contents of E are destroyed. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
          eigenvectors of the matrix A, with the i-th column of Z
          holding the eigenvector associated with D(i). \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (fla_max(1,2*N-2))
          If JOBZ = 'N', WORK is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the algorithm failed to converge; i
                off-diagonal elements of E did not converge to zero. \n

 *  * */
    template <typename T>
    void stev(char *jobz, integer *n, T *d, T *e, T *z, integer *ldz, T *work, integer *info)
    {
        stev(jobz, n, d, e, z, ldz, work, info);
    }
    /** @}*/ // end of stev

    /** @defgroup stevd stevd
     * @ingroup SYM
     * @{
     */
    /*! @brief STEVD computes the eigenvalues and, optionally, the left and/or right
 eigenvectors for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     STEVD computes all eigenvalues and, optionally, eigenvectors of a
     real symmetric tridiagonal matrix. If eigenvectors are desired, it
     uses a divide and conquer algorithm.

     The divide and conquer algorithm makes very mild assumptions about
     floating point arithmetic. It will work on machines with a guard
     digit in add/subtract, or on those binary machines without guard
     digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     Cray-2. It could conceivably fail on hexadecimal or decimal machines
     without guard digits, but we know of none.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the n diagonal elements of the tridiagonal matrix
          A.
          On exit, if INFO = 0, the eigenvalues in ascending order. \n
 * @param[in,out] E
          E is REAL array, dimension (N-1) \n
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix A, stored in elements 1 to N-1 of E.
          On exit, the contents of E are destroyed. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
          eigenvectors of the matrix A, with the i-th column of Z
          holding the eigenvector associated with D(i). \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array,
                                         dimension (LWORK) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1. \n
          If JOBZ  = 'V' and N > 1 then LWORK must be at least
                         ( 1 + 4*N + N**2 ). \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal sizes of the WORK and IWORK
          arrays, returns these values as the first entries of the WORK
          and IWORK arrays, and no error message related to LWORK or
          LIWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1. \n
          If JOBZ  = 'V' and N > 1 then LIWORK must be at least 3+5*N. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal sizes of the WORK and
          IWORK arrays, returns these values as the first entries of
          the WORK and IWORK arrays, and no error message related to
          LWORK or LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the algorithm failed to converge; i
                off-diagonal elements of E did not converge to zero. \n

 *  * */
    template <typename T>
    void stevd(char *jobz, integer *n, T *d, T *e, T *z, integer *ldz, T *work, integer *lwork,
               integer *iwork, integer *liwork, integer *info)
    {
        stevd(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
    }
    /** @}*/ // end of stevd

    /** @defgroup stevr stevr
     * @ingroup SYM
     * @{
     */
    /*! @brief STEVR computes the eigenvalues and, optionally, the left and/or right
 eigenvectors for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     STEVR computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric tridiagonal matrix T.  Eigenvalues and
     eigenvectors can be selected by specifying either a range of values
     or a range of indices for the desired eigenvalues.

     Whenever possible, SSTEVR calls SSTEMR to compute the
     eigenspectrum using Relatively Robust Representations.  SSTEMR
     computes eigenvalues by the dqds algorithm, while orthogonal
     eigenvectors are computed from various "good" L D L^T representations
     (also known as Relatively Robust Representations). Gram-Schmidt
     orthogonalization is avoided as far as possible. More specifically,
     the various steps of the algorithm are as follows. For the i-th
     unreduced block of T,
        (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
             is a relatively robust representation,
        (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
            relative accuracy by the dqds algorithm,
        (c) If there is a cluster of close eigenvalues, "choose" sigma_i
            close to the cluster, and go to step (a),
        (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
            compute the corresponding eigenvector by forming a
            rank-revealing twisted factorization.
     The desired accuracy of the output can be specified by the input
     parameter ABSTOL.

     For more details, see "A new O(n^2) algorithm for the symmetric
     tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
     Computer Science Division Technical Report No. UCB//CSD-97-971,
     UC Berkeley, May 1997.

     Note 1 : SSTEVR calls SSTEMR when the full spectrum is requested
     on machines which conform to the ieee-754 floating point standard.
     SSTEVR calls SSTEBZ and SSTEIN on non-ieee machines and
     when partial spectrum requests are made.

     Normal execution of SSTEMR may create NaNs and infinities and
     hence may abort due to a floating point exception in environments
     which do not handle NaNs and infinities in the ieee standard default
     manner.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and
          SSTEIN are called \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the n diagonal elements of the tridiagonal matrix
          A. \n
          On exit, D may be multiplied by a constant factor chosen
          to avoid over/underflow in computing the eigenvalues. \n
 * @param[in,out] E
          E is REAL array, dimension (fla_max(1,N-1)) \n
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix A in elements 1 to N-1 of E. \n
          On exit, E may be multiplied by a constant factor chosen
          to avoid over/underflow in computing the eigenvalues. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form. \n
 \n
          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3. \n
 \n
          If high relative accuracy is important, set ABSTOL to
          SLAMCH( 'Safe minimum').  Doing so will guarantee that
          eigenvalues are computed to high relative accuracy when
          possible in future releases.  The current code does not
          make any guarantees about high relative accuracy, but
          future releases will. See J. Barlow and J. Demmel,
          "Computing Accurate Eigensystems of Scaled Diagonally
          Dominant Matrices", LAPACK Working Note #7, for a discussion
          of which matrices define their eigenvalues to high relative
          accuracy. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first M elements contain the selected eigenvalues in
          ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out] ISUPPZ
          ISUPPZ is INTEGER array, dimension ( 2*max(1,M)) \n
          The support of the eigenvectors in Z, i.e., the indices
          indicating the nonzero elements in Z. The i-th eigenvector
          is nonzero only in elements ISUPPZ( 2*i-1) through
          ISUPPZ( 2*i).
          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal (and
          minimal) LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= 20*N. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal sizes of the WORK and IWORK
          arrays, returns these values as the first entries of the WORK
          and IWORK arrays, and no error message related to LWORK or
          LIWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the optimal (and
          minimal) LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK.  LIWORK >= 10*N. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal sizes of the WORK and
          IWORK arrays, returns these values as the first entries of
          the WORK and IWORK arrays, and no error message related to
          LWORK or LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  Internal error \n

 *  * */
    template <typename T>
    void stevr(char *jobz, char *range, integer *n, T *d, T *e, T *vl, T *vu, integer *il,
               integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz, integer *isuppz,
               T *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
    {
        stevr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
              iwork, liwork, info);
    }
    /** @}*/ // end of stevr

    /** @defgroup stevx stevx
     * @ingroup SYM
     * @{
     */
    /*! @brief STEVX computes the eigenvalues and, optionally, the left and/or right
eigenvectors for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     STEVX computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric tridiagonal matrix A.  Eigenvalues and
     eigenvectors can be selected by specifying either a range of values
     or a range of indices for the desired eigenvalues.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the n diagonal elements of the tridiagonal matrix
          A. \n
          On exit, D may be multiplied by a constant factor chosen
          to avoid over/underflow in computing the eigenvalues. \n
 * @param[in,out] E
          E is REAL array, dimension (fla_max(1,N-1)) \n
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix A in elements 1 to N-1 of E. \n
          On exit, E may be multiplied by a constant factor chosen
          to avoid over/underflow in computing the eigenvalues. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less
          than or equal to zero, then  EPS*|T|  will be used in
          its place, where |T| is the 1-norm of the tridiagonal
          matrix. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*SLAMCH('S'). \n
 \n
          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first M elements contain the selected eigenvalues in
          ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          If an eigenvector fails to converge (INFO > 0), then that
          column of Z contains the latest approximation to the
          eigenvector, and the index of the eigenvector is returned
          in IFAIL. \n
          If JOBZ = 'N', then Z is not referenced.
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
* @param[out]	WORK
          WORK is REAL array, dimension (5*N) \n
* @param[out]	IWORK
          IWORK is INTEGER array, dimension (5*N) \n
* @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (N) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
* @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, then i eigenvectors failed to converge.
                Their indices are stored in array IFAIL. \n

 *  * */
    template <typename T>
    void stevx(char *jobz, char *range, integer *n, T *d, T *e, T *vl, T *vu, integer *il,
               integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz, float *work,
               integer *iwork, integer *ifail, integer *info)
    {
        stevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
    }
    /** @}*/ // end of stevx

    /** @defgroup pteqr pteqr
     * @ingroup SYM
     * @{
     */
    /*! @brief PTEQR computes all eigenvalues and, optionally, eigenvectors of a symmetric
 matrix

 * @details
 * \b Purpose:
    \verbatim
     PTEQR computes all eigenvalues and, optionally, eigenvectors of a
     symmetric positive definite tridiagonal matrix by first factoring the
     matrix using SPTTRF, and then calling SBDSQR to compute the singular
     values of the bidiagonal factor.

     This routine computes the eigenvalues of the positive definite
     tridiagonal matrix to high relative accuracy.  This means that if the
     eigenvalues range over many orders of magnitude in size, then the
     small eigenvalues and corresponding eigenvectors will be computed
     more accurately than, for example, with the standard QR method.

     The eigenvectors of a full or band symmetric positive definite matrix
     can also be found if SSYTRD, SSPTRD, or SSBTRD has been used to
     reduce this matrix to tridiagonal form. (The reduction to tridiagonal
     form, however, may preclude the possibility of obtaining high
     relative accuracy in the small eigenvalues of the original matrix, if
     these eigenvalues range over many orders of magnitude.)
    \endverbatim

 * @param[in] COMPZ
          COMPZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only. \n
          = 'V':  Compute eigenvectors of original symmetric
                  matrix also.  Array Z contains the orthogonal
                  matrix used to reduce the original matrix to
                  tridiagonal form. \n
          = 'I':  Compute eigenvectors of tridiagonal matrix also. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the n diagonal elements of the tridiagonal
          matrix. \n
          On normal exit, D contains the eigenvalues, in descending
          order. \n
 * @param[in,out] E
          E is REAL array, dimension (N-1) \n
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix.
          On exit, E has been destroyed. \n
 * @param[in,out] Z
          Z is REAL array, dimension (LDZ, N) \n
          On entry, if COMPZ = 'V', the orthogonal matrix used in the
          reduction to tridiagonal form. \n
          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the
          original symmetric matrix; \n
          if COMPZ = 'I', the orthonormal eigenvectors of the
          tridiagonal matrix. \n
          If INFO > 0 on exit, Z contains the eigenvectors associated
          with only the stored eigenvalues. \n
          If  COMPZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          COMPZ = 'V' or 'I', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (4*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, and i is: \n
                <= N  the Cholesky factorization of the matrix could
                      not be performed because the i-th principal minor
                      was not positive definite. \n
                > N   the SVD algorithm failed to converge;
                      if INFO = N+i, i off-diagonal elements of the
                      bidiagonal factor did not converge to zero. \n

 *  * */
    template <typename T>
    void pteqr(char *compz, integer *n, T *d, T *e, T *z, integer *ldz, T *work, integer *info)
    {
        pteqr(compz, n, d, e, z, ldz, work, info);
    }
    template <typename T, typename Ta>
    void pteqr(char *compz, integer *n, Ta *d, Ta *e, T *z, integer *ldz, Ta *work, integer *info)
    {
        pteqr(compz, n, d, e, z, ldz, work, info);
    }
    /** @}*/ // end of pteqr

    /** @defgroup stebz stebz
     * @ingroup SYM
     * @{
     */
    /*! @brief STEBZ computes the eigenvalues of a symmetric tridiagonal matrix T

 * @details
 * \b Purpose:
    \verbatim
     STEBZ computes the eigenvalues of a symmetric tridiagonal
     matrix T.  The user may ask for all eigenvalues, all eigenvalues
     in the half-open interval (VL, VU], or the IL-th through IU-th
     eigenvalues.

     To avoid overflow, the matrix must be scaled so that its
     largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value,
 and for greatest accuracy, it should not be much smaller than that.

     See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     Matrix", Report CS41, Computer Science Dept., Stanford
     University, July 21, 1966.
    \endverbatim

 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': ("All")   all eigenvalues will be found. \n
          = 'V': ("Value") all eigenvalues in the half-open interval
                           (VL, VU] will be found. \n
          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
                           entire matrix) will be found. \n
 * @param[in] ORDER
          ORDER is CHARACTER*1 \n
          = 'B': ("By Block") the eigenvalues will be grouped by
                              split-off block (see IBLOCK, ISPLIT) and
                              ordered from smallest to largest within
                              the block. \n
          = 'E': ("Entire matrix")
                              the eigenvalues for the entire matrix
                              will be ordered from smallest to
                              largest. \n
 * @param[in] N
          N is INTEGER \n
          The order of the tridiagonal matrix T.  N >= 0. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues.  Eigenvalues less than or equal
          to VL, or greater than VU, will not be returned.  VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues.  Eigenvalues less than or equal
          to VL, or greater than VU, will not be returned.  VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'.
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute tolerance for the eigenvalues.  An eigenvalue
          (or cluster) is considered to be located if it has been
          determined to lie in an interval whose width is ABSTOL or
          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
          will be used, where |T| means the 1-norm of T. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*SLAMCH('S'), not zero. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the tridiagonal matrix T. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) off-diagonal elements of the tridiagonal matrix T. \n
 * @param[out] M
          M is INTEGER \n
          The actual number of eigenvalues found. 0 <= M <= N.
          (See also the description of INFO=2,3.) \n
 * @param[out] NSPLIT
          NSPLIT is INTEGER \n
          The number of diagonal blocks in the matrix T.
          1 <= NSPLIT <= N. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          On exit, the first M elements of W will contain the
          eigenvalues.  (SSTEBZ may use the remaining N-M elements as
          workspace.) \n
 * @param[out] IBLOCK
          IBLOCK is INTEGER array, dimension (N) \n
          At each row/column j where E(j) is zero or small, the
          matrix T is considered to split into a block diagonal
          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
          block (from 1 to the number of blocks) the eigenvalue W(i)
          belongs.  (SSTEBZ may use the remaining N-M elements as
          workspace.) \n
 * @param[out] ISPLIT
          ISPLIT is INTEGER array, dimension (N) \n
          The splitting points, at which T breaks up into submatrices.
          The first submatrix consists of rows/columns 1 to ISPLIT(1),
          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
          etc., and the NSPLIT-th consists of rows/columns
          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
          (Only the first NSPLIT elements will actually be used, but
          since the user cannot know a priori what value NSPLIT will
          have, N words must be reserved for ISPLIT.) \n
 * @param[out]	WORK
          WORK is REAL array, dimension (4*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (3*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  some or all of the eigenvalues failed to converge or
                were not computed: \n
                =1 or 3: Bisection failed to converge for some
                        eigenvalues; these eigenvalues are flagged by a
                        negative block number.  The effect is that the
                        eigenvalues may not be as accurate as the
                        absolute and relative tolerances.  This is
                        generally caused by unexpectedly inaccurate
                        arithmetic. \n
                =2 or 3: RANGE='I' only: Not all of the eigenvalues
                        IL:IU were found. \n
                        Effect: M < IU+1-IL \n
                        Cause:  non-monotonic arithmetic, causing the
                                Sturm sequence to be non-monotonic. \n
                        Cure:   recalculate, using RANGE='A', and pick
                                out eigenvalues IL:IU.  In some cases,
                                increasing the PARAMETER "FUDGE" may
                                make things work. \n
                = 4:    RANGE='I', and the Gershgorin interval
                        initially used was too small.  No eigenvalues
                        were computed. \n
                        Probable cause: your machine has sloppy
                                        floating-point arithmetic. \n
                        Cure: Increase the PARAMETER "FUDGE",
                              recompile, and try again. \n

 *  * */
    template <typename T>
    void stebz(char *range, char *order, integer *n, T *vl, T *vu, integer *il, integer *iu,
               T *abstol, T *d, T *e, integer *m, integer *nsplit, T *w, integer *iblock,
               integer *isplit, T *work, integer *iwork, integer *info)
    {
        stebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work,
              iwork, info);
    }
    /** @}*/ // end of stebz

    /** @defgroup stefr stefr
     * @ingroup SYM
     * @{
     */
    /*! @brief STERF computes all eigenvalues of a symmetric tridiagonal matrix

 * @details
 * \b Purpose:
    \verbatim
     STERF computes all eigenvalues of a symmetric tridiagonal matrix
     using the Pal-Walker-Kahan variant of the QL or QR algorithm.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the n diagonal elements of the tridiagonal matrix.
          On exit, if INFO = 0, the eigenvalues in ascending order. \n
 * @param[in,out] E
          E is REAL array, dimension (N-1) \n
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix.
          On exit, E has been destroyed. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  the algorithm failed to find all of the eigenvalues in
                a total of 30*N iterations; if INFO = i, then i
                elements of E have not converged to zero. \n

 *  * */
    template <typename T>
    void sterf(integer *n, T *d, T *e, integer *info)
    {
        sterf(n, d, e, info);
    }
    /** @}*/ // end of stefr

    /** @defgroup stedc stedc
     * @ingroup SYM
     * @{
     */
    /*! @brief Tridiagonal divide-and-conquer algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal
        matrix using the divide and conquer method. The eigenvectors of a full or band real
        symmetric matrix can also be found if SSYTRD or SSPTRD or SSBTRD has been used to reduce
        this matrix to tridiagonal form.

        This code makes very mild assumptions about floating point arithmetic. It will work on
        machines with a guard digit in add/subtract, or on those binary machines without guard
        digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could
        conceivably fail on hexadecimal or decimal machines without guard digits, but we know of
        none.  See SLAED3 for details.
    \endverbatim

    * @param[in] compz
              compz is char* \n
              = 'N':  Compute eigenvalues only. \n
              = 'I':  Compute eigenvectors of tridiagonal matrix also. \n
              = 'V':  Compute eigenvectors of original dense symmetric
              matrix also.  On entry, Z contains the orthogonal
              matrix used to reduce the original matrix to
              tridiagonal form. \n
    * @param[in] n
              n is integer* \n
              The dimension of the symmetric tridiagonal matrix.  n >= 0. \n
    * @param[in,out] d
              d is float/double array, dimension (n) \n
              On entry, the diagonal elements of the tridiagonal matrix. \n
              On exit, if info = 0, the eigenvalues in ascending order. \n
    * @param[in,out] e
              e is float/double array, dimension (n-1) \n
              On entry, the subdiagonal elements of the tridiagonal matrix. \n
              On exit, E has been destroyed. \n
    * @param[in,out] z
              z is float/double array, dimension (ldz,n) \n
              On entry, if compz = 'V', then Z contains the orthogonal
              matrix used in the reduction to tridiagonal form. \n
              On exit, if info = 0, then if compz = 'V', Z contains the
              orthonormal eigenvectors of the original symmetric matrix,
              and if compz = 'I', Z contains the orthonormal eigenvectors
              of the symmetric tridiagonal matrix. \n
              If  compz = 'N', then Z is not referenced. \n
    * @param[in] ldz
              ldz is integer* \n
              The leading dimension of the array Z.  ldz >= 1. \n
              If eigenvectors are desired, then ldz >= fla_max(1,n). \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1. \n
              If COMPZ = 'V' and N > 1, LWORK must be at least N*N. \n
              Note that for COMPZ = 'V', then if N is less than or
              equal to the minimum divide size, usually 25, then LWORK need
              only be 1. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal sizes of the WORK, RWORK and
              IWORK arrays, returns these values as the first entries of
              the WORK, RWORK and IWORK arrays, and no error message
              related to LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	RWORK
              RWORK is REAL array, dimension (MAX(1,LRWORK)) \n
              On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. \n
    * @param[in]	LRWORK
              LRWORK is INTEGER \n
              The dimension of the array RWORK. \n
              If COMPZ = 'N' or N <= 1, LRWORK must be at least 1. \n
              If COMPZ = 'V' and N > 1, LRWORK must be at least \n
                             1 + 3*N + 2*N*lg N + 4*N**2 , \n
                             where lg( N ) = smallest integer k such \n
                             that 2**k >= N. \n
              If COMPZ = 'I' and N > 1, LRWORK must be at least
                             1 + 4*N + 2*N**2 . \n
              Note that for COMPZ = 'I' or 'V', then if N is less than or
              equal to the minimum divide size, usually 25, then LRWORK
              need only be fla_max(1,2*(N-1)). \n
 \n
              If LRWORK = -1, then a workspace query is assumed; the
              routine only calculates the optimal sizes of the WORK, RWORK
              and IWORK arrays, returns these values as the first entries
              of the WORK, RWORK and IWORK arrays, and no error message
              related to LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	IWORK
              IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
              On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
    * @param[in]	LIWORK
              LIWORK is INTEGER \n
              The dimension of the array IWORK. \n
              If COMPZ = 'N' or N <= 1, LIWORK must be at least 1. \n
              If COMPZ = 'V' or N > 1,  LIWORK must be at least
                                        6 + 6*N + 5*N*lg N. \n
              If COMPZ = 'I' or N > 1,  LIWORK must be at least
                                        3 + 5*N . \n
              Note that for COMPZ = 'I' or 'V', then if N is less than or
              equal to the minimum divide size, usually 25, then LIWORK
              need only be 1. \n
 \n
              If LIWORK = -1, then a workspace query is assumed; the
              routine only calculates the optimal sizes of the WORK, RWORK
              and IWORK arrays, returns these values as the first entries
              of the WORK, RWORK and IWORK arrays, and no error message
              related to LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit. \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n
              > 0:  The algorithm failed to compute an eigenvalue while
                    working on the submatrix lying in rows and columns
                    INFO/(N+1) through mod(INFO,N+1). \n

    *     *  */
    template <typename T>
    void stedc(char *compz, integer *n, T *d, T *e, T *z, integer *ldz, T *work, integer *lwork,
               integer *iwork, integer *liwork, integer *info)
    {
        stedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
    }
    template <typename T, typename Ta>
    void stedc(char *compz, integer *n, Ta *d, Ta *e, T *z, integer *ldz, T *work, integer *lwork,
               Ta *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info)
    {
        stedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
    }
    /** @}*/ // end of stedc

    /** @defgroup stegr stegr
     * @ingroup SYM
     * @{
     */
    /*! @brief STEGR computes selected eigenvalues and, optionally, eigenvectors of a real
 symmetric tridiagonal matrix T

 * @details
 * \b Purpose:
    \verbatim
     STEGR computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     a well defined set of pairwise different real eigenvalues, the corresponding
     real eigenvectors are pairwise orthogonal.

     The spectrum may be computed either completely or partially by specifying
     either an interval (VL,VU] or a range of indices IL:IU for the desired
     eigenvalues.

     SSTEGR is a compatibility wrapper around the improved SSTEMR routine.
     See SSTEMR for further details.

     One important change is that the ABSTOL parameter no longer provides any
     benefit and hence is no longer used.

     Note : SSTEGR and SSTEMR work only on machines which follow
     IEEE-754 floating-point standard in their handling of infinities and
     NaNs.  Normal execution may create these exceptiona values and hence
     may abort due to a floating point exception in environments which
     do not conform to the IEEE-754 standard.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the N diagonal elements of the tridiagonal matrix
          T. On exit, D is overwritten. \n
 * @param[in,out] E
          E is REAL array, dimension (N) \n
          On entry, the (N-1) subdiagonal elements of the tridiagonal
          matrix T in elements 1 to N-1 of E. E(N) need not be set on
          input, but is used internally as workspace. \n
          On exit, E is overwritten. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          Unused.  Was the absolute error tolerance for the
          eigenvalues/eigenvectors in previous versions. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first M elements contain the selected eigenvalues in
          ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z
          contain the orthonormal eigenvectors of the matrix T
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i). \n
          If JOBZ = 'N', then Z is not referenced. \n
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used.
          Supplying N columns is always safe. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', then LDZ >= fla_max(1,N). \n
 * @param[out] ISUPPZ
          ISUPPZ is INTEGER array, dimension ( 2*max(1,M)) \n
          The support of the eigenvectors in Z, i.e., the indices
          indicating the nonzero elements in Z. The i-th computed eigenvector
          is nonzero only in elements ISUPPZ( 2*i-1) through
          ISUPPZ( 2*i). This is relevant in the case when the matrix
          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (LWORK) \n
          On exit, if INFO = 0, WORK(1) returns the optimal
          (and minimal) LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,18*N)
          if JOBZ = 'V', and LWORK >= fla_max(1,12*N) if JOBZ = 'N'. \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (LIWORK) \n
          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK.  LIWORK >= fla_max(1,10*N)
          if the eigenvectors are desired, and LIWORK >= fla_max(1,8*N)
          if only the eigenvalues are to be computed. \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the IWORK array,
          returns this value as the first entry of the IWORK array, and
          no error message related to LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          On exit, INFO \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = 1X, internal error in SLARRE,
                if INFO = 2X, internal error in SLARRV.
                Here, the digit X = ABS( IINFO ) < 10, where IINFO is
                the nonzero error code returned by SLARRE or
                SLARRV, respectively. \n

 *  * */
    template <typename T>
    void stegr(char *jobz, char *range, integer *n, T *d, T *e, T *vl, T *vu, integer *il,
               integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz, integer *isuppz,
               T *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
    {
        stegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
              iwork, liwork, info);
    }
    template <typename T, typename Ta>
    void stegr(char *jobz, char *range, integer *n, Ta *d, Ta *e, Ta *vl, Ta *vu, integer *il,
               integer *iu, Ta *abstol, integer *m, Ta *w, T *z, integer *ldz, integer *isuppz,
               Ta *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
    {
        stegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
              iwork, liwork, info);
    }
    /** @}*/ // end of stegr

    /** @defgroup stein stein
     * @ingroup SYM
     * @{
     */
    /*! @brief STEIN computes the eigenvectors of a real symmetric tridiagonal matrix T

 * @details
 * \b Purpose:
    \verbatim
     STEIN computes the eigenvectors of a real symmetric tridiagonal
     matrix T corresponding to specified eigenvalues, using inverse
     iteration.

     The maximum number of iterations allowed for each eigenvector is
     specified by an internal parameter MAXITS (currently set to 5).
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the tridiagonal matrix T. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of the tridiagonal matrix
          T, in elements 1 to N-1. \n
 * @param[in] M
          M is INTEGER \n
          The number of eigenvectors to be found.  0 <= M <= N. \n
 * @param[in] W
          W is REAL array, dimension (N) \n
          The first M elements of W contain the eigenvalues for
          which eigenvectors are to be computed.  The eigenvalues
          should be grouped by split-off block and ordered from
          smallest to largest within the block.  ( The output array
          W from SSTEBZ with ORDER = 'B' is expected here.) \n
 * @param[in] IBLOCK
          IBLOCK is INTEGER array, dimension (N) \n
          The submatrix indices associated with the corresponding
          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
          the first submatrix from the top, =2 if W(i) belongs to
          the second submatrix, etc.  ( The output array IBLOCK
          from SSTEBZ is expected here.) \n
 * @param[in] ISPLIT
          ISPLIT is INTEGER array, dimension (N) \n
          The splitting points, at which T breaks up into submatrices.
          The first submatrix consists of rows/columns 1 to
          ISPLIT( 1), the second of rows/columns ISPLIT( 1)+1
          through ISPLIT( 2), etc.
          ( The output array ISPLIT from SSTEBZ is expected here.) \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, M) \n
          The computed eigenvectors.  The eigenvector associated
          with the eigenvalue W(i) is stored in the i-th column of
          Z.  Any vector which fails to converge is set to its current
          iterate after MAXITS iterations. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (5*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (M) \n
          On normal exit, all elements of IFAIL are zero.
          If one or more eigenvectors fail to converge after
          MAXITS iterations, then their indices are stored in
          array IFAIL. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit. \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, then i eigenvectors failed to converge
               in MAXITS iterations.  Their indices are stored in
               array IFAIL. \n

 *  * */
    template <typename T>
    void stein(integer *n, T *d, T *e, integer *m, T *w, integer *iblock, integer *isplit, T *z,
               integer *ldz, T *work, integer *iwork, integer *ifail, integer *info)
    {
        stein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    }
    template <typename T, typename Ta>
    void stein(integer *n, Ta *d, Ta *e, integer *m, Ta *w, integer *iblock, integer *isplit, T *z,
               integer *ldz, Ta *work, integer *iwork, integer *ifail, integer *info)
    {
        stein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    }
    /** @}*/ // end of stein

    /** @defgroup stemr stemr
     * @ingroup SYM
     * @{
     */
    /*! @brief STEMR computes selected eigenvalues and, optionally, eigenvectors of a real symmetric tridiagonal matrix T.

 * @details
 * \b Purpose:
    \verbatim
      STEMR computes selected eigenvalues and, optionally, eigenvectors
      of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
      a well defined set of pairwise different real eigenvalues, the corresponding
      real eigenvectors are pairwise orthogonal.
      The spectrum may be computed either completely or partially by specifying
      either an interval (VL,VU] or a range of indices IL:IU for the desired
      eigenvalues.
      Depending on the number of desired eigenvalues, these are computed either
      by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
      computed by the use of various suitable L D L^T factorizations near clusters
      of close eigenvalues (referred to as RRRs, Relatively Robust
      Representations). An informal sketch of the algorithm follows.
      For each unreduced block (submatrix) of T,
         (a) Compute T - sigma I  = L D L^T, so that L and D
            define all the wanted eigenvalues to high relative accuracy.
            This means that small relative changes in the entries of D and L
            cause only small relative changes in the eigenvalues and
            eigenvectors. The standard (unfactored) representation of the
            tridiagonal matrix T does not have this property in general.
         (b) Compute the eigenvalues to suitable accuracy.
            If the eigenvectors are desired, the algorithm attains full
            accuracy of the computed eigenvalues only right before
            the corresponding vectors have to be computed, see steps c) and d).
         (c) For each cluster of close eigenvalues, select a new
            shift close to the cluster, find a new factorization, and refine
            the shifted eigenvalues to suitable accuracy.
         (d) For each eigenvalue with a large enough relative separation compute
            the corresponding eigenvector by forming a rank revealing twisted
            factorization. Go back to (c) for any clusters that remain.
    \endverbatim

 * @param[in] JOBZ
          JOBS is CHARACTER*1 \n
          = 'N': Compute eigenvalues only \n
          = 'V': Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': All eigenvalues will be found. \n
          = 'V': All eigenvalues in the half-open interval (VL,VU] will be found. \n
          = 'I': The IL-th through IU-th eigenvalues will be found. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, The n diagonal elements of the tridiagonal matrix T. \n
          On exit, D is overwritten \n
 * @param[in,out] E
          E is REAL array, dimension (N) \n
          On entry, The (n-1) subdiagonal elements of the tridiagonal matrix 
          T, in elements 1 to N-1 of E. E(N) need not to be set on input,
          but is used internally as workspace. \n
          On exit, E is overwritten. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to be searched for 
          eigenvalues. VL < VU \n
          Not referenced if RANGE ='A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to be searched for 
          eigenvalues. VL < VU \n
          Not referenced if RANGE ='A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0. \n
          Not referenced if RANGE ='A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0. \n
          Not referenced if RANGE ='A' or 'V'. \n
 * @param[out] M
          M is INTEGER \n
          The number of eigenvectors to be found.  0 <= M <= N. \n
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1 \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first M elements contain the selected eigenvalues in 
          ascending order. \n
 * @param[out] Z
         Z is COMPLEX array, dimension (LDZ, max(1,M) ) \n
         If JOBZ = 'V', and if INFO = 0, then the first M columns of Z
         contain the orthonormal eigenvectors of the matrix T
         corresponding to the selected eigenvalues, with the i-th
         column of Z holding the eigenvector associated with W(i). \n
         If JOBZ = 'N', then Z is not referenced. \n
         Note: the user must ensure that at least max(1,M) columns are
         supplied in the array Z; if RANGE = 'V', the exact value of M
         is not known in advance and can be computed with a workspace
         query by setting NZC = -1, see below. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if 
          JOBZ = 'V', then LDZ >= max(1,N). \n
 * @param[in] NZC
         NZC is INTEGER \n
         The number of eigenvectors to be held in the array Z. \n
         If RANGE = 'A', then NZC >= max(1,N). \n
         If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU]. \n
         If RANGE = 'I', then NZC >= IU-IL+1. \n
         If NZC = -1, then a workspace query is assumed; the
         routine calculates the number of columns of the array Z that
         are needed to hold the eigenvectors. \n
         This value is returned as the first entry of the Z array, and
         no error message related to NZC is issued by XERBLA. \n
 * @param[out] ISUPPZ
         ISUPPZ is INTEGER array, dimension ( 2*max(1,M) ) \n
         The support of the eigenvectors in Z, i.e., the indices
         indicating the nonzero elements in Z. The i-th computed eigenvector
         is nonzero only in elements ISUPPZ( 2*i-1 ) through
         ISUPPZ( 2*i ). This is relevant in the case when the matrix
         is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0. \n
 * @param[in,out] TRYRAC
         TRYRAC is LOGICAL \n
         If TRYRAC = .TRUE., indicates that the code should check whether
         the tridiagonal matrix defines its eigenvalues to high relative
         accuracy.  If so, the code uses relative-accuracy preserving
         algorithms that might be (a bit) slower depending on the matrix. \n
         If the matrix does not define its eigenvalues to high relative
         accuracy, the code can uses possibly faster algorithms. \n
         If TRYRAC = .FALSE., the code is not required to guarantee
         relatively accurate eigenvalues and can use the fastest possible
         techniques. \n
         On exit, a .TRUE. TRYRAC will be set to .FALSE. if the matrix
         does not define its eigenvalues to high relative accuracy. \n
 * @param[out] WORK
         WORK is REAL array, dimension (LWORK) \n
         On exit, if INFO = 0, WORK(1) returns the optimal
         (and minimal) LWORK. \n
 * @param[in] LWORK
         LWORK is INTEGER \n
         The dimension of the array WORK. LWORK >= max(1,18*N) \n
         if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'. \n
         If LWORK = -1, then a workspace query is assumed; the routine
         only calculates the optimal size of the WORK array, returns
         this value as the first entry of the WORK array, and no error
         message related to LWORK is issued by XERBLA. \n
 * @param[out] IWORK
         IWORK is INTEGER array, dimension (LIWORK) \n
         On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
 * @param[in] LIWORK
         LIWORK is INTEGER \n
         The dimension of the array IWORK.  LIWORK >= max(1,10*N) \n
         if the eigenvectors are desired, and LIWORK >= max(1,8*N)
         if only the eigenvalues are to be computed. \n
         If LIWORK = -1, then a workspace query is assumed; the
         routine only calculates the optimal size of the IWORK array,
         returns this value as the first entry of the IWORK array, and
         no error message related to LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
         INFO is INTEGER \n
         On exit, INFO \n
         = 0:  successful exit \n
         < 0:  if INFO = -i, the i-th argument had an illegal value \n
         > 0:  if INFO = 1X, internal error in SLARRE, \n
         if INFO = 2X, internal error in CLARRV. \n
         Here, the digit X = ABS( IINFO ) < 10, where IINFO is
         the nonzero error code returned by SLARRE or
         CLARRV, respectively. \n

 *  * */
   template <typename T>
   void stemr(char* jobz, char* range, integer* n, T*  d, T*  e, T* vl, T* vu, integer* il,
              integer* iu, integer* m, T*  w, T* z, integer* ldz, integer* nzc, integer* isuppz,
              integer* tryrac, T* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
   {
      stemr(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac,
            work, lwork, iwork, liwork, info);
   }
   template< typename T, typename Ta >
   void stemr(char* jobz, char* range, integer* n, Ta*  d, Ta*  e, Ta* vl, Ta* vu, integer* il,
              integer* iu, integer* m, Ta*  w, T* z, integer* ldz, integer* nzc, integer* isuppz,
              integer* tryrac, Ta* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
   {
      stemr(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac,
            work, lwork, iwork, liwork, info);
   }
    /** @}*/ // end of stemr

    /** @defgroup steqr steqr
     * @ingroup SYM
     * @{
     */
    /*! @brief Tridiagonal QR algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal
        matrix using the implicit QL or QR method. The eigenvectors of a full or band symmetric
        matrix can also be found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this
    matrix to tridiagonal form. \endverbatim

    * @param[in] jobz
              jobz is char* \n
              = 'N':  Compute eigenvalues only. \n
              = 'V':  Compute eigenvalues and eigenvectors of the original
              symmetric matrix.  On entry, Z must contain the
              orthogonal matrix used to reduce the original matrix
              to tridiagonal form. \n
              = 'I':  Compute eigenvalues and eigenvectors of the
              tridiagonal matrix.  Z is initialized to the identity
              matrix. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix.  n >= 0. \n
    * @param[in,out] d
              d is float/double array, dimension (n) \n
              On entry, the diagonal elements of the tridiagonal matrix. \n
              On exit, if info = 0, the eigenvalues in ascending order. \n
    * @param[in,out] e
              e is float/double array, dimension (n-1) \n
              On entry, the (n-1) subdiagonal elements of the tridiagonal
              matrix. \n
              On exit, E has been destroyed. \n
    * @param[in,out] z
              z is float/double array, dimension (ldz, n) \n
              On entry, if  jobz = 'V', then Z contains the orthogonal
              matrix used in the reduction to tridiagonal form. \n
              On exit, if info = 0, then if  jobz = 'V', Z contains the
              orthonormal eigenvectors of the original symmetric matrix,
              and if jobz = 'I', Z contains the orthonormal eigenvectors
              of the symmetric tridiagonal matrix. \n
              If jobz = 'N', then Z is not referenced. \n
    * @param[in] ldz
              ldz is integer* \n
              The leading dimension of the array Z.  ldz >= 1, and if
              eigenvectors are desired, then  ldz >= fla_max(1,n). \n
    * @param[out]	WORK
              WORK is REAL array, dimension (fla_max(1,2*N-2)) \n
              If COMPZ = 'N', then WORK is not referenced. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  the algorithm has failed to find all the eigenvalues in
                    a total of 30*N iterations; if INFO = i, then i
                    elements of E have not converged to zero; on exit, D
                    and E contain the elements of a symmetric tridiagonal
                    matrix which is unitarily similar to the original
                    matrix. \n

    *     *  */
    template <typename T>
    void steqr(char *compz, integer *n, T *d, T *e, T *z, integer *ldz, T *work, integer *info)
    {
        steqr(compz, n, d, e, z, ldz, work, info);
    }
    template <typename T, typename Ta>
    void steqr(char *compz, integer *n, Ta *d, Ta *e, T *z, integer *ldz, Ta *work, integer *info)
    {
        steqr(compz, n, d, e, z, ldz, work, info);
    }
    /** @}*/ // end of steqr

    /** @defgroup sygv sygv
     * @ingroup SYM
     * @{
     */
    /*! @brief SYGV  computes all the eigenvalues, the eigenvectors \n
     of a real generalized symmetric-definite eigenproblem
 * @details
 * \b Purpose:
    \verbatim
     SYGV computes all the eigenvalues, and optionally, the eigenvectors
     of a real generalized symmetric-definite eigenproblem, of the form
     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
     Here A and B are assumed to be symmetric and B is also
     positive definite.
    \endverbatim

 * @param[in] ITYPE
          ITYPE is INTEGER \n
          Specifies the problem type to be solved: \n
          = 1:  A*x = (lambda)*B*x \n
          = 2:  A*B*x = (lambda)*x \n
          = 3:  B*A*x = (lambda)*x \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A. \n
 \n
          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
          matrix Z of eigenvectors.  The eigenvectors are normalized
          as follows: \n
          if ITYPE = 1 or 2, Z**T*B*Z = I; \n
          if ITYPE = 3, Z**T*inv(B)*Z = I. \n
          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
          or the lower triangle (if UPLO='L') of A, including the
          diagonal, is destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the symmetric positive definite matrix B. \n
          If UPLO = 'U', the leading N-by-N upper triangular part of B
          contains the upper triangular part of the matrix B. \n
          If UPLO = 'L', the leading N-by-N lower triangular part of B
          contains the lower triangular part of the matrix B. \n
 \n
          On exit, if INFO <= N, the part of B containing the matrix is
          overwritten by the triangular factor U or L from the Cholesky
          factorization B = U**T*U or B = L*L**T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of the array WORK.  LWORK >= fla_max(1,3*N-1).
          For optimal efficiency, LWORK >= (NB+2)*N,
          where NB is the blocksize for SSYTRD returned by ILAENV. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  SPOTRF or SSYEV returned an error code: \n
             <= N:  if INFO = i, SSYEV failed to converge;
                    i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero; \n
             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                    minor of order i of B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void sygv(integer *itype, char *jobz, char *uplo, integer *n, T *a, integer *lda, T *b,
              integer *ldb, T *w, T *work, integer *lwork, integer *info)
    {
        sygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
    }
    /** @}*/ // end of sygv

    /** @defgroup sygv_2stage sygv_2stage
     * @ingroup SYM
     * @{
     */
    /*! @brief SYGV_2STAGE  computes all the eigenvalues, the eigenvectors \n
     of a real generalized symmetric-definite eigenproblem
 * @details
 * \b Purpose:
    \verbatim
     SYGV_2STAGE computes all the eigenvalues, and optionally, the eigenvectors
     of a real generalized symmetric-definite eigenproblem, of the form
     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
     Here A and B are assumed to be symmetric and B is also
     positive definite.
     This routine use the 2stage technique for the reduction to tridiagonal
     which showed higher performance on recent architecture and for large
     sizes N>2000.
    \endverbatim

 * @param[in] ITYPE
          ITYPE is INTEGER \n
          Specifies the problem type to be solved: \n
          = 1:  A*x = (lambda)*B*x \n
          = 2:  A*B*x = (lambda)*x \n
          = 3:  B*A*x = (lambda)*x \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors.
                  Not available in this release. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A. \n
 \n
          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
          matrix Z of eigenvectors.  The eigenvectors are normalized
          as follows: \n
          if ITYPE = 1 or 2, Z**T*B*Z = I; \n
          if ITYPE = 3, Z**T*inv(B)*Z = I. \n
          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
          or the lower triangle (if UPLO='L') of A, including the
          diagonal, is destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the symmetric positive definite matrix B.
          If UPLO = 'U', the leading N-by-N upper triangular part of B
          contains the upper triangular part of the matrix B.
          If UPLO = 'L', the leading N-by-N lower triangular part of B
          contains the lower triangular part of the matrix B. \n
 \n
          On exit, if INFO <= N, the part of B containing the matrix is
          overwritten by the triangular factor U or L from the Cholesky
          factorization B = U**T*U or B = L*L**T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of the array WORK. LWORK >= 1, when N <= 1; \n
          otherwise   \n
          If JOBZ = 'N' and N > 1, LWORK must be queried.
                                   LWORK = MAX(1, dimension) where
                                   dimension = fla_max(stage1,stage2) + (KD+1)*N + 2*N
                                             = N*KD + N*max(KD+1,FACTOPTNB)
                                               + fla_max(2*KD*KD, KD*NTHREADS)
                                               + (KD+1)*N + 2*N
                                   where KD is the blocking size of the reduction,
                                   FACTOPTNB is the blocking used by the QR or LQ
                                   algorithm, usually FACTOPTNB=128 is a good choice
                                   NTHREADS is the number of threads used when
                                   openMP compilation is enabled, otherwise =1. \n
          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  SPOTRF or SSYEV returned an error code: \n
             <= N:  if INFO = i, SSYEV failed to converge;
                    i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero; \n
             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                    minor of order i of B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void sygv_2stage(integer *itype, char *jobz, char *uplo, integer *n, T *a, integer *lda, T *b,
                     integer *ldb, T *w, T *work, integer *lwork, integer *info)
    {
        sygv_2stage(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
    }
    /** @}*/ // end of sygv_2stage

    /** @defgroup sygvd sygvd
     * @ingroup SYM
     * @{
     */
    /*! @brief SYGVD  computes all the eigenvalues, the eigenvectors \n
     of a real generalized symmetric-definite eigenproblem
 * @details
 * \b Purpose:
    \verbatim
     SYGVD computes all the eigenvalues, and optionally, the eigenvectors
     of a real generalized symmetric-definite eigenproblem, of the form
     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
     B are assumed to be symmetric and B is also positive definite.
     If eigenvectors are desired, it uses a divide and conquer algorithm.

     The divide and conquer algorithm makes very mild assumptions about
     floating point arithmetic. It will work on machines with a guard
     digit in add/subtract, or on those binary machines without guard
     digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     Cray-2. It could conceivably fail on hexadecimal or decimal machines
     without guard digits, but we know of none.
    \endverbatim

 * @param[in] ITYPE
          ITYPE is INTEGER \n
          Specifies the problem type to be solved: \n
          = 1:  A*x = (lambda)*B*x \n
          = 2:  A*B*x = (lambda)*x \n
          = 3:  B*A*x = (lambda)*x \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A. \n
 \n
          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
          matrix Z of eigenvectors.  The eigenvectors are normalized
          as follows: \n
          if ITYPE = 1 or 2, Z**T*B*Z = I; \n
          if ITYPE = 3, Z**T*inv(B)*Z = I. \n
          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
          or the lower triangle (if UPLO='L') of A, including the
          diagonal, is destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the symmetric matrix B.  If UPLO = 'U', the
          leading N-by-N upper triangular part of B contains the
          upper triangular part of the matrix B.  If UPLO = 'L',
          the leading N-by-N lower triangular part of B contains
          the lower triangular part of the matrix B. \n
 \n
          On exit, if INFO <= N, the part of B containing the matrix is
          overwritten by the triangular factor U or L from the Cholesky
          factorization B = U**T*U or B = L*L**T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order.
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If N <= 1,               LWORK >= 1. \n
          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1. \n
          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal sizes of the WORK and IWORK
          arrays, returns these values as the first entries of the WORK
          and IWORK arrays, and no error message related to LWORK or
          LIWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          If N <= 1,                LIWORK >= 1. \n
          If JOBZ  = 'N' and N > 1, LIWORK >= 1. \n
          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal sizes of the WORK and
          IWORK arrays, returns these values as the first entries of
          the WORK and IWORK arrays, and no error message related to
          LWORK or LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  SPOTRF or SSYEVD returned an error code: \n
             <= N:  if INFO = i and JOBZ = 'N', then the algorithm
                    failed to converge; i off-diagonal elements of an
                    intermediate tridiagonal form did not converge to
                    zero;
                    if INFO = i and JOBZ = 'V', then the algorithm
                    failed to compute an eigenvalue while working on
                    the submatrix lying in rows and columns INFO/(N+1)
                    through mod(INFO,N+1); \n
             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                    minor of order i of B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void sygvd(integer *itype, char *jobz, char *uplo, integer *n, T *a, integer *lda, T *b,
               integer *ldb, T *w, T *work, integer *lwork, integer *iwork, integer *liwork,
               integer *info)
    {
        sygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info);
    }
    /** @}*/ // end of sygvd

    /** @defgroup sygvx sygvx
     * @ingroup SYM
     * @{
     */
    /*! @brief SYGVX computes selected eigenvalues, the eigenvectors \n
     of a real generalized symmetric-definite eigenproblem
 * @details
 * \b Purpose:
    \verbatim
     SYGVX computes selected eigenvalues, and optionally, eigenvectors
     of a real generalized symmetric-definite eigenproblem, of the form
     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
     and B are assumed to be symmetric and B is also positive definite.
     Eigenvalues and eigenvectors can be selected by specifying either a
     range of values or a range of indices for the desired eigenvalues.
    \endverbatim

 * @param[in] ITYPE
          ITYPE is INTEGER \n
          Specifies the problem type to be solved: \n
          = 1:  A*x = (lambda)*B*x \n
          = 2:  A*B*x = (lambda)*x \n
          = 3:  B*A*x = (lambda)*x \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A and B are stored; \n
          = 'L':  Lower triangle of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix pencil (A,B).  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA, N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A. \n
 \n
          On exit, the lower triangle (if UPLO='L') or the upper
          triangle (if UPLO='U') of A, including the diagonal, is
          destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, N) \n
          On entry, the symmetric matrix B.  If UPLO = 'U', the
          leading N-by-N upper triangular part of B contains the
          upper triangular part of the matrix B.  If UPLO = 'L',
          the leading N-by-N lower triangular part of B contains
          the lower triangular part of the matrix B. \n
 \n
          On exit, if INFO <= N, the part of B containing the matrix is
          overwritten by the triangular factor U or L from the Cholesky
          factorization B = U**T*U or B = L*L**T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing C to tridiagonal form, where C is the symmetric
          matrix of the standard symmetric problem to which the
          generalized problem is transformed. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*SLAMCH('S'). \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          On normal exit, the first M elements contain the selected
          eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'N', then Z is not referenced. \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i). \n
          The eigenvectors are normalized as follows: \n
          if ITYPE = 1 or 2, Z**T*B*Z = I; \n
          if ITYPE = 3, Z**T*inv(B)*Z = I. \n
 \n
          If an eigenvector fails to converge, then that column of Z
          contains the latest approximation to the eigenvector, and the
          index of the eigenvector is returned in IFAIL. \n
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out] IFAIL
          IFAIL is INTEGER array, dimension (N) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of the array WORK.  LWORK >= fla_max(1,8*N). \n
          For optimal efficiency, LWORK >= (NB+3)*N,
          where NB is the blocksize for SSYTRD returned by ILAENV. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (5*N) \n
 * @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (N) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  SPOTRF or SSYEVX returned an error code: \n
             <= N:  if INFO = i, SSYEVX failed to converge;
                    i eigenvectors failed to converge.  Their indices
                    are stored in array IFAIL. \n
             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                    minor of order i of B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void sygvx(integer *itype, char *jobz, char *range, char *uplo, integer *n, T *a, integer *lda,
               T *b, integer *ldb, T *vl, T *vu, integer *il, integer *iu, T *abstol, integer *m,
               T *w, T *z, integer *ldz, T *work, integer *lwork, integer *iwork, integer *ifail,
               integer *info)
    {
        sygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz,
              work, lwork, iwork, ifail, info);
    }
    /** @}*/ // end of sygvx

    /** @defgroup spgv spgv
     * @ingroup SYM
     * @{
     */
    /*! @brief SPGV computes all the eigenvalues

 * @details
 * \b Purpose:
    \verbatim
     SPGV computes all the eigenvalues and, optionally, the eigenvectors
     of a real generalized symmetric-definite eigenproblem, of the form
     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
     Here A and B are assumed to be symmetric, stored in packed format,
     and B is also positive definite.
    \endverbatim

 * @param[in] ITYPE
          ITYPE is INTEGER \n
          Specifies the problem type to be solved: \n
          = 1:  A*x = (lambda)*B*x \n
          = 2:  A*B*x = (lambda)*x \n
          = 3:  B*A*x = (lambda)*x \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, the contents of AP are destroyed. \n
 * @param[in,out] BP
          BP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          B, packed columnwise in a linear array.  The j-th column of B
          is stored in the array BP as follows: \n
          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; \n
          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. \n
 \n
          On exit, the triangular factor U or L from the Cholesky
          factorization B = U**T*U or B = L*L**T, in the same storage
          format as B. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
          eigenvectors.  The eigenvectors are normalized as follows:
          if ITYPE = 1 or 2, Z**T*B*Z = I; \n
          if ITYPE = 3, Z**T*inv(B)*Z = I. \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  SPPTRF or SSPEV returned an error code: \n
             <= N:  if INFO = i, SSPEV failed to converge;
                    i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero. \n
             > N:   if INFO = n + i, for 1 <= i <= n, then the leading
                    minor of order i of B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void spgv(integer *itype, char *jobz, char *uplo, integer *n, T *ap, T *bp, T *w, T *z,
              integer *ldz, T *work, integer *info)
    {
        spgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info);
    }
    template <typename T, typename Ta>
    void hpgv(integer *itype, char *jobz, char *uplo, integer *n, T *ap, T *bp, T *w, T *z,
              integer *ldz, T *work, Ta *rwork, integer *info)
    {
        hpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, info);
    }
    /** @}*/ // end of spgv

    /** @defgroup spgvd spgvd
     * @ingroup SYM
     * @{
     */
    /*! @brief SPGVD computes all the eigenvalues

 * @details
 * \b Purpose:
    \verbatim
     SPGVD computes all the eigenvalues, and optionally, the eigenvectors
     of a real generalized symmetric-definite eigenproblem, of the form
     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
     B are assumed to be symmetric, stored in packed format, and B is also
     positive definite.
     If eigenvectors are desired, it uses a divide and conquer algorithm.

     The divide and conquer algorithm makes very mild assumptions about
     floating point arithmetic. It will work on machines with a guard
     digit in add/subtract, or on those binary machines without guard
     digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     Cray-2. It could conceivably fail on hexadecimal or decimal machines
     without guard digits, but we know of none.
    \endverbatim

 * @param[in] ITYPE \n
          ITYPE is INTEGER \n
          Specifies the problem type to be solved: \n
          = 1:  A*x = (lambda)*B*x \n
          = 2:  A*B*x = (lambda)*x \n
          = 3:  B*A*x = (lambda)*x \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, the contents of AP are destroyed. \n
 * @param[in,out] BP
          BP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          B, packed columnwise in a linear array.  The j-th column of B
          is stored in the array BP as follows: \n
          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; \n
          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. \n
 \n
          On exit, the triangular factor U or L from the Cholesky
          factorization B = U**T*U or B = L*L**T, in the same storage
          format as B. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
          eigenvectors.  The eigenvectors are normalized as follows: \n
          if ITYPE = 1 or 2, Z**T*B*Z = I; \n
          if ITYPE = 3, Z**T*inv(B)*Z = I. \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the required LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If N <= 1,               LWORK >= 1. \n
          If JOBZ = 'N' and N > 1, LWORK >= 2*N. \n
          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the required sizes of the WORK and IWORK
          arrays, returns these values as the first entries of the WORK
          and IWORK arrays, and no error message related to LWORK or
          LIWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if INFO = 0, IWORK(1) returns the required LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          If JOBZ  = 'N' or N <= 1, LIWORK >= 1. \n
          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the required sizes of the WORK and
          IWORK arrays, returns these values as the first entries of
          the WORK and IWORK arrays, and no error message related to
          LWORK or LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  SPPTRF or SSPEVD returned an error code: \n
             <= N:  if INFO = i, SSPEVD failed to converge;
                    i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero; \n
             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                    minor of order i of B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void spgvd(integer *itype, char *jobz, char *uplo, integer *n, T *ap, T *bp, T *w, T *z,
               integer *ldz, T *work, integer *lwork, integer *iwork, integer *liwork,
               integer *info)
    {
        spgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, iwork, liwork, info);
    }
    template <typename T, typename Ta>
    void hpgvd(integer *itype, char *jobz, char *uplo, integer *n, T *ap, T *bp, Ta *w, T *z,
               integer *ldz, T *work, integer *lwork, Ta *rwork, integer *lrwork, integer *iwork,
               integer *liwork, integer *info)
    {
        hpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork,
              info);
    }
    /** @}*/ // end of spgvd

    /** @defgroup spgvx spgvx
     * @ingroup SYM
     * @{
     */
    /*! @brief SPGVX computes all the eigenvalues

 * @details
 * \b Purpose:
    \verbatim
     SPGVX computes selected eigenvalues, and optionally, eigenvectors
     of a real generalized symmetric-definite eigenproblem, of the form
     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
     and B are assumed to be symmetric, stored in packed storage, and B
     is also positive definite.  Eigenvalues and eigenvectors can be
     selected by specifying either a range of values or a range of indices
     for the desired eigenvalues.
    \endverbatim

 * @param[in] ITYPE
          ITYPE is INTEGER \n
          Specifies the problem type to be solved: \n
          = 1:  A*x = (lambda)*B*x \n
          = 2:  A*B*x = (lambda)*x \n
          = 3:  B*A*x = (lambda)*x \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A and B are stored; \n
          = 'L':  Lower triangle of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix pencil (A,B).  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, the contents of AP are destroyed. \n
 * @param[in,out] BP
          BP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          B, packed columnwise in a linear array.  The j-th column of B
          is stored in the array BP as follows: \n
          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; \n
          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. \n
 \n
          On exit, the triangular factor U or L from the Cholesky
          factorization B = U**T*U or B = L*L**T, in the same storage
          format as B. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*SLAMCH('S'). \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          On normal exit, the first M elements contain the selected
          eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If JOBZ = 'N', then Z is not referenced. \n
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          The eigenvectors are normalized as follows: \n
          if ITYPE = 1 or 2, Z**T*B*Z = I; \n
          if ITYPE = 3, Z**T*inv(B)*Z = I. \n
 \n
          If an eigenvector fails to converge, then that column of Z
          contains the latest approximation to the eigenvector, and the
          index of the eigenvector is returned in IFAIL.
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (8*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (5*N) \n
 * @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (N) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  SPPTRF or SSPEVX returned an error code: \n
             <= N:  if INFO = i, SSPEVX failed to converge;
                    i eigenvectors failed to converge.  Their indices
                    are stored in array IFAIL. \n
             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                    minor of order i of B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void spgvx(integer *itype, char *jobz, char *range, char *uplo, integer *n, T *ap, T *bp, T *vl,
               T *vu, integer *il, integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz,
               T *work, integer *iwork, integer *ifail, integer *info)
    {
        spgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work,
              iwork, ifail, info);
    }
    template <typename T, typename Ta>
    void hpgvx(integer *itype, char *jobz, char *range, char *uplo, integer *n, T *ap, T *bp,
               Ta *vl, Ta *vu, integer *il, integer *iu, Ta *abstol, integer *m, Ta *w, T *z,
               integer *ldz, T *work, Ta *rwork, integer *iwork, integer *ifail, integer *info)
    {
        hpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work,
              rwork, iwork, ifail, info);
    }
    /** @}*/ // end of spgvx

    /** @defgroup sbgv sbgv
     * @ingroup SYM
     * @{
     */
    /*! @brief SBGV computes all the eigenvalues, and optionally, the eigenvectors

 * @details
 * \b Purpose:
    \verbatim
     SBGV computes all the eigenvalues, and optionally, the eigenvectors
     of a real generalized symmetric-definite banded eigenproblem, of
     the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric
     and banded, and B is also positive definite.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in] KA
          KA is INTEGER \n
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'. KA >= 0. \n
 * @param[in] KB
          KB is INTEGER \n
          The number of superdiagonals of the matrix B if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'. KB >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB, N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix A, stored in the first ka+1 rows of the array.  The
          j-th column of A is stored in the j-th column of the array AB
          as follows: \n
          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for fla_max(1,j-ka)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). \n
 \n
          On exit, the contents of AB are destroyed. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KA+1. \n
 * @param[in,out] BB
          BB is REAL array, dimension (LDBB, N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix B, stored in the first kb+1 rows of the array.  The
          j-th column of B is stored in the j-th column of the array BB
          as follows: \n
          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for fla_max(1,j-kb)<=i<=j; \n
          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). \n
 \n
          On exit, the factor S from the split Cholesky factorization
          B = S**T*S, as returned by SPBSTF. \n
 * @param[in] LDBB
          LDBB is INTEGER \n
          The leading dimension of the array BB.  LDBB >= KB+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
          eigenvectors, with the i-th column of Z holding the
          eigenvector associated with W(i). The eigenvectors are
          normalized so that Z**T*B*Z = I. \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= N. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, and i is: \n
             <= N:  the algorithm failed to converge:
                    i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero; \n
             > N:   if INFO = N + i, for 1 <= i <= N, then SPBSTF
                    returned INFO = i: B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void sbgv(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, T *ab, integer *ldab,
              T *bb, integer *ldbb, T *w, T *z, integer *ldz, float *work, integer *info)
    {
        sbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, info);
    }
    template <typename T, typename Ta>
    void hbgv(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, T *ab, integer *ldab,
              T *bb, integer *ldbb, Ta *w, T *z, integer *ldz, T *work, Ta *rwork, integer *info)
    {
        hbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info);
    }
    /** @}*/ // end of sbgv

    /** @defgroup sbgvd sbgvd
     * @ingroup SYM
     * @{
     */
    /*! @brief SBGVD computes all the eigenvalues, and optionally, the eigenvectors

 * @details
 * \b Purpose:
    \verbatim
     SBGVD computes all the eigenvalues, and optionally, the eigenvectors
     of a real generalized symmetric-definite banded eigenproblem, of the
     form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and
     banded, and B is also positive definite.  If eigenvectors are
     desired, it uses a divide and conquer algorithm.

     The divide and conquer algorithm makes very mild assumptions about
     floating point arithmetic. It will work on machines with a guard
     digit in add/subtract, or on those binary machines without guard
     digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     Cray-2. It could conceivably fail on hexadecimal or decimal machines
     without guard digits, but we know of none.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in] KA
          KA is INTEGER \n
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KA >= 0. \n
 * @param[in] KB
          KB is INTEGER \n
          The number of superdiagonals of the matrix B if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KB >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB, N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix A, stored in the first ka+1 rows of the array.  The
          j-th column of A is stored in the j-th column of the array AB
          as follows: \n
          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for fla_max(1,j-ka)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). \n
 \n
          On exit, the contents of AB are destroyed. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KA+1. \n
 * @param[in,out] BB
          BB is REAL array, dimension (LDBB, N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix B, stored in the first kb+1 rows of the array.  The
          j-th column of B is stored in the j-th column of the array BB
          as follows: \n
          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for fla_max(1,j-kb)<=i<=j; \n
          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). \n
 \n
          On exit, the factor S from the split Cholesky factorization
          B = S**T*S, as returned by SPBSTF. \n
 * @param[in] LDBB
          LDBB is INTEGER \n
          The leading dimension of the array BB.  LDBB >= KB+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
          eigenvectors, with the i-th column of Z holding the
          eigenvector associated with W(i).  The eigenvectors are
          normalized so Z**T*B*Z = I. \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If N <= 1,               LWORK >= 1. \n
          If JOBZ = 'N' and N > 1, LWORK >= 3*N. \n
          If JOBZ = 'V' and N > 1, LWORK >= 1 + 5*N + 2*N**2. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal sizes of the WORK and IWORK
          arrays, returns these values as the first entries of the WORK
          and IWORK arrays, and no error message related to LWORK or
          LIWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          If JOBZ  = 'N' or N <= 1, LIWORK >= 1. \n
          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal sizes of the WORK and
          IWORK arrays, returns these values as the first entries of
          the WORK and IWORK arrays, and no error message related to
          LWORK or LIWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, and i is: \n
             <= N:  the algorithm failed to converge:
                    i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero; \n
             > N:   if INFO = N + i, for 1 <= i <= N, then SPBSTF
                    returned INFO = i: B is not positive definite.
                    The factorization of B could not be completed and
                    no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void sbgvd(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, T *ab, integer *ldab,
               T *bb, integer *ldbb, T *w, T *z, integer *ldz, T *work, integer *lwork,
               integer *iwork, integer *liwork, integer *info)
    {
        sbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork, liwork,
              info);
    }
    template <typename T, typename Ta>
    void hbgvd(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, T *ab, integer *ldab,
               T *bb, integer *ldbb, Ta *w, T *z, integer *ldz, T *work, integer *lwork, Ta *rwork,
               integer *lrwork, integer *iwork, integer *liwork, integer *info)
    {
        hbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, rwork, lrwork,
              iwork, liwork, info);
    }
    /** @}*/ // end of sbgvd

    /** @defgroup sgbvx sgbvx
     * @ingroup SYM
     * @{
     */
    /*! @brief SBGVX computes all the eigenvalues, and optionally, the eigenvectors

 * @details
 * \b Purpose:
    \verbatim
     SBGVX computes selected eigenvalues, and optionally, eigenvectors
     of a real generalized symmetric-definite banded eigenproblem, of
     the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric
     and banded, and B is also positive definite.  Eigenvalues and
     eigenvectors can be selected by specifying either all eigenvalues,
     a range of values or a range of indices for the desired eigenvalues.
    \endverbatim

 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute eigenvalues only; \n
          = 'V':  Compute eigenvalues and eigenvectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all eigenvalues will be found. \n
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th eigenvalues will be found. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangles of A and B are stored; \n
          = 'L':  Lower triangles of A and B are stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in] KA
          KA is INTEGER \n
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KA >= 0. \n
 * @param[in] KB
          KB is INTEGER \n
          The number of superdiagonals of the matrix B if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KB >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB, N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix A, stored in the first ka+1 rows of the array.  The
          j-th column of A is stored in the j-th column of the array AB
          as follows: \n
          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for fla_max(1,j-ka)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). \n
 \n
          On exit, the contents of AB are destroyed. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KA+1. \n
 * @param[in,out] BB
          BB is REAL array, dimension (LDBB, N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix B, stored in the first kb+1 rows of the array.  The
          j-th column of B is stored in the j-th column of the array BB
          as follows: \n
          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for fla_max(1,j-kb)<=i<=j; \n
          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). \n
 \n
          On exit, the factor S from the split Cholesky factorization
          B = S**T*S, as returned by SPBSTF. \n
 * @param[in] LDBB
          LDBB is INTEGER \n
          The leading dimension of the array BB.  LDBB >= KB+1. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ, N) \n
          If JOBZ = 'V', the n-by-n matrix used in the reduction of
          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,
          and consequently C to tridiagonal form. \n
          If JOBZ = 'N', the array Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  If JOBZ = 'N',
          LDQ >= 1. If JOBZ = 'V', LDQ >= fla_max(1,N). \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to \n
 \n
                  ABSTOL + EPS *   fla_max( |a|,|b|) , \n
 \n
          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form. \n
 \n
          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*SLAMCH('S'). \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues found.  0 <= M <= N. \n
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          If INFO = 0, the eigenvalues in ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, N) \n
          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
          eigenvectors, with the i-th column of Z holding the
          eigenvector associated with W(i).  The eigenvectors are
          normalized so Z**T*B*Z = I. \n
          If JOBZ = 'N', then Z is not referenced. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (7*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (5*N) \n
 * @param[out]	IFAIL
          IFAIL is INTEGER array, dimension (M) \n
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvalues that failed to converge. \n
          If JOBZ = 'N', then IFAIL is not referenced. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          <= N: if INFO = i, then i eigenvectors failed to converge.
                  Their indices are stored in IFAIL. \n
          > N:  SPBSTF returned an error code; i.e.,
                if INFO = N + i, for 1 <= i <= N, then the leading
                minor of order i of B is not positive definite.
                The factorization of B could not be completed and
                no eigenvalues or eigenvectors were computed. \n

 *  * */
    template <typename T>
    void sbgvx(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb, T *ab,
               integer *ldab, T *bb, integer *ldbb, T *q, integer *ldq, T *vl, T *vu, integer *il,
               integer *iu, T *abstol, integer *m, T *w, T *z, integer *ldz, T *work,
               integer *iwork, integer *ifail, integer *info)
    {
        sbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m,
              w, z, ldz, work, iwork, ifail, info);
    }
    template <typename T, typename Ta>
    void hbgvx(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb, T *ab,
               integer *ldab, T *bb, integer *ldbb, T *q, integer *ldq, Ta *vl, Ta *vu, integer *il,
               integer *iu, Ta *abstol, integer *m, Ta *w, T *z, integer *ldz, T *work, Ta *rwork,
               integer *iwork, integer *ifail, integer *info)
    {
        hbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m,
              w, z, ldz, work, rwork, iwork, ifail, info);
    }
    /** @}*/ // end of sgbvx

    /** @defgroup disna disna
     * @ingroup SYM
     * @{
     */
    /*! @brief DISNA computes the reciprocal condition numbers for the eigenvectors \n
     of a real symmetric or complex Hermitian matrix
 * @details
 * \b Purpose:
    \verbatim
     DISNA computes the reciprocal condition numbers for the eigenvectors
     of a real symmetric or complex Hermitian matrix or for the left or
     right singular vectors of a general m-by-n matrix. The reciprocal
     condition number is the 'gap' between the corresponding eigenvalue or
     singular value and the nearest other one.

     The bound on the error, measured by angle in radians, in the I-th
     computed vector is given by

            SLAMCH( 'E') * ( ANORM / SEP( I))

     where ANORM = 2-norm(A) = fla_max( abs( D(j))).  SEP(I) is not allowed
     to be smaller than SLAMCH( 'E')*ANORM in order to limit the size of
     the error bound.

     SDISNA may also be used to compute error bounds for eigenvectors of
     the generalized symmetric definite eigenproblem.
    \endverbatim

 * @param[in] JOB
          JOB is CHARACTER*1 \n
          Specifies for which problem the reciprocal condition numbers
          should be computed: \n
          = 'E':  the eigenvectors of a symmetric/Hermitian matrix; \n
          = 'L':  the left singular vectors of a general matrix; \n
          = 'R':  the right singular vectors of a general matrix. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          If JOB = 'L' or 'R', the number of columns of the matrix,
          in which case N >= 0. Ignored if JOB = 'E'. \n
 * @param[in] D
          D is REAL array, dimension (M) if JOB = 'E'
                              dimension (min(M,N)) if JOB = 'L' or 'R' \n
          The eigenvalues (if JOB = 'E') or singular values (if JOB =
          'L' or 'R') of the matrix, in either increasing or decreasing
          order. If singular values, they must be non-negative. \n
 * @param[out] SEP
          SEP is REAL array, dimension (M) if JOB = 'E'
                               dimension (min(M,N)) if JOB = 'L' or 'R' \n
          The reciprocal condition numbers of the vectors. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void disna(char *job, integer *m, integer *n, T *d, T *sep, integer *info)
    {
        disna(job, m, n, d, sep, info);
    }
    /** @}*/ // end of disna

    /** @defgroup sytrd sytrd
     * @ingroup SYM
     * @{
     */
    /*! @brief Reduction of a real symmetric matrix a to real symmetric tridiagonal form
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric matrix a to real symmetric tridiagonal form T by an
        orthogonal similarity transformation:
            Q**T * A * Q = T.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the symmetric matrix a. \n
              If uplo = 'U', the leading
              n-by-n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced. \n
              If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if uplo = 'U', the diagonal and first superdiagonal
              of A are overwritten by the corresponding elements of the
              tridiagonal matrix T, and the elements above the first
              superdiagonal, with the array tau, represent the orthogonal
              matrix Q as a product of elementary reflectors; \n
              if uplo = 'L', the diagonal and first subdiagonal of A are over-
              written by the corresponding elements of the tridiagonal
              matrix T, and the elements below the first subdiagonal, with
              the array tau, represent the orthogonal matrix Q as a product
              of elementary reflectors. See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] d
              d is float/double array, dimension (n) \n
              The diagonal elements of the tridiagonal matrix T:
              D(i) = A(i,i). \n
    * @param[out] e
              e is float/double array, dimension (n-1) \n
              The off-diagonal elements of the tridiagonal matrix T: \n
              E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'. \n
    * @param[out] tau
              tau is float/double array, dimension (n-1) \n
              The scalar factors of the elementary reflectors (see Further
              Details).
              *
              * \n
              * **Further Details**
              * \verbatim
                      If uplo = 'U', the matrix Q is represented as a product of elementary
    reflectors

                      Q = H(n-1) . . . H(2) H(1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**T

                      where tau is a real scalar, and V is a real vector with V(i+1:n) = 0 and
    V(i) = 1; V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

                      If uplo = 'L', the matrix Q is represented as a product of elementary
    reflectors

                      Q = H(1) H(2) . . . H(n-1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**T

                      where tau is a real scalar, and V is a real vector with V(1:i) = 0 and
    V(i+1) = 1; V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

                      The contents of A on exit are illustrated by the following examples with n
    = 5:

                      if uplo = 'U':  if uplo = 'L':

                      (  d   e   v2  v3  v4)  (  d  )
                      (   d   e   v3  v4)  (  e   d )
                      ( d   e   v4)  (  v1  e   d)
                      (  d   e )  (  v1  v2  e   d  )
                      (   d )  (  v1  v2  v3  e   d )

                      where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i). \endverbatim
    * @param[out]	WORK
              WORK is REAL array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK.  LWORK >= 1. \n
              For optimum performance LWORK >= N*NB, where NB is the
              optimal blocksize. \n

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
    void sytrd(char *uplo, integer *n, T *a, integer *lda, T *d, T *e, T *tau, T *work,
               integer *lwork, integer *info)
    {
        sytrd(uplo, n, a, lda, d, e, tau, work, lwork, info);
    }
    /** @}*/ // end of sytrd

    /** @defgroup hetrd hetrd
     * @ingroup SYM
     * @{
     */
    /*! @brief Reduction of a complex Hermitian matrix a to real symmetric tridiagonal form
    *
    * @details
    * \b Purpose:
    * \verbatim
    Reduction of a complex Hermitian matrix a to real symmetric tridiagonal form T by a
    unitary similarity transformation:
        Q**H * A * Q = T.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the symmetric matrix a. \n
              If uplo = 'U', the leading
              n-by-n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced. \n
              If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if uplo = 'U', the diagonal and first superdiagonal
              of A are overwritten by the corresponding elements of the
              tridiagonal matrix T, and the elements above the first
              superdiagonal, with the array tau, represent the orthogonal
              matrix Q as a product of elementary reflectors; \n
              if uplo = 'L', the diagonal and first subdiagonal of A are over-
              written by the corresponding elements of the tridiagonal
              matrix T, and the elements below the first subdiagonal, with
              the array tau, represent the orthogonal matrix Q as a product
              of elementary reflectors. See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] d
              d is float/double array, dimension (n) \n
              The diagonal elements of the tridiagonal matrix T:
              D(i) = A(i,i). \n
    * @param[out] e
              e is float/double array, dimension (n-1) \n
              The off-diagonal elements of the tridiagonal matrix T: \n
              E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'. \n
    * @param[out] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (n-1) \n
              The scalar factors of the elementary reflectors (see Further
              Details).
              *
              * \n
              * **Further Details**
              * \verbatim
                      If uplo = 'U', the matrix Q is represented as a product of elementary
    reflectors

                      Q = H(n-1) . . . H(2) H(1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**H

                      where tau is a complex scalar, and V is a complex vector with V(i+1:n) = 0
    and V(i) = 1; V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

                      If uplo = 'L', the matrix Q is represented as a product of elementary
    reflectors

                      Q = H(1) H(2) . . . H(n-1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**H

                      where tau is a complex scalar, and V is a complex vector with V(1:i) = 0
    and V(i+1) = 1; V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

                      The contents of A on exit are illustrated by the following examples with n
    = 5:

                      if uplo = 'U':  if uplo = 'L':

                      (  d   e   v2  v3  v4)  (  d  )
                      (   d   e   v3  v4)  (  e   d )
                      ( d   e   v4)  (  v1  e   d)
                      (  d   e )  (  v1  v2  e   d  )
                      (   d )  (  v1  v2  v3  e   d )

                      where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i). \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK.  LWORK >= 1.
              For optimum performance LWORK >= N*NB, where NB is the
              optimal blocksize. \n

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T, typename Ta>
    void hetrd(char *uplo, integer *n, T *a, integer *lda, Ta *d, Ta *e, T *tau, T *work,
               integer *lwork, integer *info)
    {
        hetrd(uplo, n, a, lda, d, e, tau, work, lwork, info);
    }
    /** @}*/ // end of hetrd

    /** @defgroup hetd2 hetd2
     * @ingroup SYM
     * @{
     */
    /*! @brief Reduction of a Hermitian matrix a to real symmetric tridiagonal form (unblocked
    algorithm)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a Hermitian matrix a to real symmetric tridiagonal form T by an unitary
        similarity transformation(unblocked algorithm):
            Q**T * A * Q = T.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the symmetric matrix a. \n
              If uplo = 'U', the leading
              n-by-n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced. \n
              If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if uplo = 'U', the diagonal and first superdiagonal
              of A are overwritten by the corresponding elements of the
              tridiagonal matrix T, and the elements above the first
              superdiagonal, with the array tau, represent the orthogonal
              matrix Q as a product of elementary reflectors; \n
              if uplo = 'L', the diagonal and first subdiagonal of A are over-
              written by the corresponding elements of the tridiagonal
              matrix T, and the elements below the first subdiagonal, with
              the array tau, represent the orthogonal matrix Q as a product
              of elementary reflectors. See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] d
              d is float/double array, dimension (n) \n
              The diagonal elements of the tridiagonal matrix T:
              D(i) = A(i,i).
    * @param[out] e
              e is float/double array, dimension (n-1) \n
              The off-diagonal elements of the tridiagonal matrix T: \n
              E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'. \n
    * @param[out] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (n-1) \n
              The scalar factors of the elementary reflectors (see Further
              Details).
              *
              * \n
              * **Further Details**
              * \verbatim
                      If uplo = 'U', the matrix Q is represented as a product of elementary
    reflectors

                          Q = H(n-1) . . . H(2) H(1).

                      Each H(i) has the form

                          H(i) = I - tau * V * V**H

                      where tau is a complex scalar, and V is a complex vector with V(i+1:n) = 0
    and V(i) = 1; V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

                      If uplo = 'L', the matrix Q is represented as a product of elementary
    reflectors

                          Q = H(1) H(2) . . . H(n-1).

                      Each H(i) has the form

                          H(i) = I - tau * V * V**H

                      where tau is a complex scalar, and V is a complex complex vector with
    V(1:i) = 0 and V(i+1) = 1; V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

                      The contents of A on exit are illustrated by the following examples with n
    = 5:

                      if uplo = 'U':                       if uplo = 'L':

                          (  d   e   v2  v3  v4)              (  d                 )
                          (      d   e   v3  v4)              (  e   d             )
                          (          d   e   v4)              (  v1  e   d         )
                          (              d   e )              (  v1  v2  e   d     )
                          (                  d )              (  v1  v2  v3  e   d )

                      where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i). \endverbatim
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T, typename Ta>
    void hetd2(char *uplo, integer *n, T *a, integer *lda, Ta *d, Ta *e, T *tau, integer *info)
    {
        hetd2(uplo, n, a, lda, d, e, tau, info);
    }
    /** @}*/ // end of hetd2

    /** @defgroup sytd2 sytd2
     * @ingroup SYM
     * @{
     */
    /*! @brief Reduction of a real symmetric matrix a to real symmetric tridiagonal form
    (unblocked algorithm)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric matrix a to real symmetric tridiagonal form T by an
    orthogonal similarity transformation(unblocked algorithm): Q**T * A * Q = T. \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the symmetric matrix a. \n
              If uplo = 'U', the leading
              n-by-n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced. \n
              If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if uplo = 'U', the diagonal and first superdiagonal
              of A are overwritten by the corresponding elements of the
              tridiagonal matrix T, and the elements above the first
              superdiagonal, with the array tau, represent the orthogonal
              matrix Q as a product of elementary reflectors; \n
              if uplo = 'L', the diagonal and first subdiagonal of A are over-
              written by the corresponding elements of the tridiagonal
              matrix T, and the elements below the first subdiagonal, with
              the array tau, represent the orthogonal matrix Q as a product
              of elementary reflectors. See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[out] d
              d is float/double array, dimension (n) \n
              The diagonal elements of the tridiagonal matrix T:
              D(i) = A(i,i). \n
    * @param[out] e
              e is float/double array, dimension (n-1) \n
              The off-diagonal elements of the tridiagonal matrix T: \n
              E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'. \n
    * @param[out] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (n-1) \n
              The scalar factors of the elementary reflectors (see Further
              Details).
              *
              * \n
              * **Further Details**
              * \verbatim
                      If uplo = 'U', the matrix Q is represented as a product of elementary
    reflectors

                      Q = H(n-1) . . . H(2) H(1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**T

                      where tau is a real scalar, and V is a real vector with V(i+1:n) = 0 and
    V(i) = 1; V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

                      If uplo = 'L', the matrix Q is represented as a product of elementary
    reflectors

                      Q = H(1) H(2) . . . H(n-1).

                      Each H(i) has the form

                      H(i) = I - tau * V * V**T

                      where tau is a real scalar, and V is a real vector with V(1:i) = 0 and
    V(i+1) = 1; V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

                      The contents of A on exit are illustrated by the following examples with n
    = 5:

                      if uplo = 'U':  if uplo = 'L':

                      (  d   e   v2  v3  v4)  (  d  )
                      (   d   e   v3  v4)  (  e   d )
                      ( d   e   v4)  (  v1  e   d)
                      (  d   e )  (  v1  v2  e   d  )
                      (   d )  (  v1  v2  v3  e   d )

                      where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i). \endverbatim
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void sytd2(char *uplo, integer *n, T *a, integer *lda, T *d, T *e, T *tau, integer *info)
    {
        sytd2(uplo, n, a, lda, d, e, tau, info);
    }

    /** @}*/ // end of sytd2

    /** @defgroup latrd latrd
     * @ingroup SYM
     * @{
     */
    /*! @brief LATRD reduces the first nb rows and columns of a symmetric/Hermitian   \n
     matrix A to real tridiagonal form by an orthogonal similarity transformation
 * @details
 * \b Purpose:
    \verbatim
    LATRD reduces NB rows and columns of a real symmetric matrix A to
    symmetric tridiagonal form by an orthogonal similarity
    transformation Q**T * A * Q, and   returns the matrices V and W which are
    needed to apply the transformation to the unreduced part of A.

    If UPLO = 'U', SLATRD reduces the last NB rows and columns of a
    matrix, of which the upper triangle is supplied;
    if UPLO = 'L', SLATRD reduces the first NB rows and columns of a
    matrix, of which the lower triangle is supplied.

    This is an auxiliary routine called by SSYTRD.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U': Upper triangular \n
          = 'L': Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. \n
 * @param[in] NB
          NB is INTEGER \n
          The number of rows and columns to be reduced. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          n-by-n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n-by-n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
          On exit: \n
          if UPLO = 'U', the last NB columns have been reduced to
            tridiagonal form, with the diagonal elements overwriting
            the diagonal elements of A; the elements above the diagonal
            with the array TAU, represent the orthogonal matrix Q as a
            product of elementary reflectors; \n
          if UPLO = 'L', the first NB columns have been reduced to
            tridiagonal form, with the diagonal elements overwriting
            the diagonal elements of A; the elements below the diagonal
            with the array TAU, represent the  orthogonal matrix Q as a
            product of elementary reflectors. \n
          See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= (1,N). \n
 * @param[out] E
          E is REAL array, dimension (N-1) \n
          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
          elements of the last NB columns of the reduced matrix; \n
          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
          the first NB columns of the reduced matrix. \n
 * @param[out] TAU
          TAU is REAL array, dimension (N-1) \n
          The scalar factors of the elementary reflectors, stored in
          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
          See Further Details. \n
 * @param[out] W
          W is REAL array, dimension (LDW,NB) \n
          The n-by-nb matrix W required to update the unreduced part
          of A. \n
 * @param[in] LDW
          LDW is INTEGER \n
          The leading dimension of the array W. LDW >= fla_max(1,N).  \n

 *  * */
    template <typename T>
    void latrd(char *uplo, integer *n, integer *nb, T *a, integer *lda, T *e, T *tau, T *w,
               integer *ldw)
    {
        latrd(uplo, n, nb, a, lda, e, tau, w, ldw);
    }
    template <typename T, typename Ta>
    void latrd(char *uplo, integer *n, integer *nb, T *a, integer *lda, Ta *e, T *tau, T *w,
               integer *ldw)
    {
        latrd(uplo, n, nb, a, lda, e, tau, w, ldw);
    }
    /** @}*/ // end of latrd

    /** @defgroup ungtr ungtr
     * @ingroup SYM
     * @{
     */
    /*! @brief Form Q from tridiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
    Form Q from tridiagonal reduction. Generate a complex unitary matrix Q which is defined as
    the product of n-1 elementary reflectors of order M, as returned by CHETRD:

        if uplo = 'U', Q = H(n-1) . . . H(2) H(1),

        if uplo = 'L', Q = H(1) H(2) . . . H(n-1).
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U': Upper triangle of a contains elementary reflectors
              from CHETRD; \n
              = 'L': Lower triangle of a contains elementary reflectors
              from CHETRD. \n
    * @param[in] m
              m is integer* \n
              The order of the matrix Q. m >= 0. \n
    * @param[in,out] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,m) \n
              On entry, the vectors which define the elementary reflectors,
              as returned by CHETRD. \n
              On exit, the m-by-m orthogonal matrix Q. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. lda >= fla_max(1,m). \n
    * @param[in] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (m-1) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by CHETRD. \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= N-1. \n
              For optimum performance LWORK >= (N-1)*NB, where NB is
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

    *     *  */
    template <typename T>
    void ungtr(char *uplo, integer *m, T *a, integer *lda, T *tau, T *work, integer *lwork,
               integer *info)
    {
        ungtr(uplo, m, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of ungtr

    /** @defgroup orgtr orgtr
     * @ingroup SYM
     * @{
     */
    /*! @brief Form Q from tridiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
    Form Q from tridiagonal reduction. Generate a real orthogonal matrix Q which is defined as
    the product of n-1 elementary reflectors of order M, as returned by SSYTRD:

        if uplo = 'U', Q = H(n-1) . . . H(2) H(1),

        if uplo = 'L', Q = H(1) H(2) . . . H(n-1).
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U': Upper triangle of a contains elementary reflectors
              from SSYTRD; \n
              = 'L': Lower triangle of a contains elementary reflectors
              from SSYTRD. \n
    * @param[in] m
              m is integer* \n
              The order of the matrix Q. m >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,m) \n
              On entry, the vectors which define the elementary reflectors,
              as returned by SSYTRD. \n
              On exit, the m-by-m orthogonal matrix Q. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. lda >= fla_max(1,m). \n
    * @param[in] tau
              tau is float/double array, dimension (m-1) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by SSYTRD. \n
    * @param[out]	WORK
              WORK is REAL array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= fla_max(1,N-1). \n
              For optimum performance LWORK >= (N-1)*NB, where NB is
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

    *     *  */
    template <typename T>
    void orgtr(char *uplo, integer *m, T *a, integer *lda, T *tau, T *work, integer *lwork,
               integer *info)
    {
        orgtr(uplo, m, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of orgtr

    /** @defgroup unmtr unmtr
     * @ingroup SYM
     * @{
     */
    /*! @brief Apply Q or Q' from tridiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from tridiagonal reduction. Overwrite the general complex m-by-n matrix c
 with

            side = 'L'  side = 'R'
            trans = 'N':   Q * C C * Q
            trans = 'C':   Q**H * C C * Q**H

        where Q is a complex unitary matrix of order nq, with nq = m if side = 'L' and nq = n if
        side = 'R'. Q is defined as the product of nq-1 elementary reflectors, as returned by
 CHETRD:

            if uplo = 'U', Q = H(nq-1) . . . H(2) H(1);

            if uplo = 'L', Q = H(1) H(2) . . . H(nq-1).
    \endverbatim

    * @param[in] side
              side is char* \n
              = 'L': apply Q or Q**H from the Left; \n
              = 'R': apply Q or Q**H from the Right. \n
    * @param[in] uplo
              uplo is char* \n
              = 'U': Upper triangle of a contains elementary reflectors
              from CHETRD; \n
              = 'L': Lower triangle of a contains elementary reflectors
              from CHETRD. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q; \n
              = 'C':  Transpose, apply Q**C. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] a
              a is COMPLEX/COMPLEX*16 array, dimension \n
              (lda,m) if side = 'L' \n
              (lda,n) if side = 'R' \n
              The vectors which define the elementary reflectors, as
              returned by CHETRD. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              lda >= fla_max(1,m) if side = 'L'; lda >= fla_max(1,n) if side = 'R'. \n
    * @param[in] tau
              tau is COMPLEX/COMPLEX*16 array, dimension \n
              (m-1) if side = 'L' \n
              (n-1) if side = 'R' \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by CHETRD. \n
    * @param[in,out] c
              c is COMPLEX/COMPLEX*16 array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. \n
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
              if SIDE = 'R', LWORK >= fla_max(1,M). \n
              For optimum performance LWORK >= N*NB if SIDE = 'L', and
              LWORK >=M*NB if SIDE = 'R', where NB is the optimal
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
    void unmtr(char *side, char *uplo, char *trans, integer *m, integer *n, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @}*/ // end of unmtr

    /** @defgroup ormtr ormtr
     * @ingroup SYM
     * @{
     */
    /*! @brief Apply Q or Q' from tridiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from tridiagonal reduction. Overwrite the general real m-by-n matrix c
 with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'T':   Q**T * C C * Q**T

        where Q is a real orthogonal matrix of order nq, with nq = m if side = 'L' and nq = n if
        side = 'R'. Q is defined as the product of nq-1 elementary reflectors, as returned by
 SSYTRD:

        if uplo = 'U', Q = H(nq-1) . . . H(2) H(1);

        if uplo = 'L', Q = H(1) H(2) . . . H(nq-1).
    \endverbatim

    * @param[in] side
              side is char* \n
              = 'L': apply Q or Q**T from the Left; \n
              = 'R': apply Q or Q**T from the Right. \n
    * @param[in] uplo
              uplo is char* \n
              = 'U': Upper triangle of a contains elementary reflectors
              from SSYTRD; \n
              = 'L': Lower triangle of a contains elementary reflectors
              from SSYTRD. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q; \n
              = 'T':  Transpose, apply Q**T. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] a
              a is float/double array, dimension \n
              (lda,m) if side = 'L' \n
              (lda,n) if side = 'R' \n
              The vectors which define the elementary reflectors, as
              returned by SSYTRD. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              lda >= fla_max(1,m) if side = 'L'; lda >= fla_max(1,n) if side = 'R'. \n
    * @param[in] tau
              tau is float/double array, dimension \n
              (m-1) if side = 'L' \n
              (n-1) if side = 'R' \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i), as returned by SSYTRD. \n
    * @param[in,out] c
              c is float/double array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
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
    void ormtr(char *side, char *uplo, char *trans, integer *m, integer *n, T *a, integer *lda,
               T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
    }

    /** @}*/ // end of- ormtr

    /** @defgroup sytrd_2stage sytrd_2stage
     * @ingroup SYM
     * @{
     */
    /*! @brief SYTRD_2STAGE reduces a real symmetric matrix A to real symmetric tridiagonal form
 T

 * @details
 * \b Purpose:
    \verbatim
    SYTRD_2STAGE reduces a real symmetric matrix A to real symmetric
    tridiagonal form T by a orthogonal similarity transformation:
    Q1**T Q2**T* A * Q2 * Q1 = T.
    \endverbatim

 * @param[in] VECT
          VECT is CHARACTER*1 \n
          = 'N':  No need for the Housholder representation,
                  in particular for the second stage (Band to
                  tridiagonal) and thus LHOUS2 is of size fla_max(1, 4*N); \n
          = 'V':  the Householder representation is needed to
                  either generate Q1 Q2 or to apply Q1 Q2,
                  then LHOUS2 is to be queried and computed.
                  (NOT AVAILABLE IN THIS RELEASE). \n
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
          triangular part of A is not referenced. \n \n
          On exit, if UPLO = 'U', the band superdiagonal
          of A are overwritten by the corresponding elements of the
          internal band-diagonal matrix AB, and the elements above
          the KD superdiagonal, with the array TAU, represent the orthogonal
          matrix Q1 as a product of elementary reflectors; if UPLO
          = 'L', the diagonal and band subdiagonal of A are over-
          written by the corresponding elements of the internal band-diagonal
          matrix AB, and the elements below the KD subdiagonal, with
          the array TAU, represent the orthogonal matrix Q1 as a product
          of elementary reflectors. See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] D
          D is REAL array, dimension (N) \n
          The diagonal elements of the tridiagonal matrix T. \n
 * @param[out] E
          E is REAL array, dimension (N-1) \n
          The off-diagonal elements of the tridiagonal matrix T. \n
 * @param[out] TAU
          TAU is REAL array, dimension (N-KD) \n
          The scalar factors of the elementary reflectors of
          the first stage (see Further Details). \n
 * @param[out] HOUS2
          HOUS2 is REAL array, dimension (LHOUS2) \n
          Stores the Householder representation of the stage2
          band to tridiagonal. \n
 * @param[in] LHOUS2
          LHOUS2 is INTEGER \n
          The dimension of the array HOUS2.
          If LWORK = -1, or LHOUS2 = -1,
          then a query is assumed; the routine
          only calculates the optimal size of the HOUS2 array, returns
          this value as the first entry of the HOUS2 array, and no error
          message related to LHOUS2 is issued by XERBLA. \n
          If VECT='N', LHOUS2 = fla_max(1, 4*n); \n
          if VECT='V', option not yet available. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK = MAX(1, dimension) \n
          If LWORK = -1, or LHOUS2=-1,
          then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
          LWORK = MAX(1, dimension) where \n
          dimension   = fla_max(stage1,stage2) + (KD+1)*N \n
                      = N*KD + N*max(KD+1,FACTOPTNB)  \n
                        + fla_max(2*KD*KD, KD*NTHREADS)  \n
                        + (KD+1)*N  \n
          where KD is the blocking size of the reduction,
          FACTOPTNB is the blocking used by the QR or LQ
          algorithm, usually FACTOPTNB=128 is a good choice
          NTHREADS is the number of threads used when
          openMP compilation is enabled, otherwise =1. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sytrd_2stage(char *vect, char *uplo, integer *n, T *a, integer *lda, T *d, T *e, T *tau,
                      T *hous2, integer *lhous2, T *work, integer *lwork, integer *info)
    {
        sytrd_2stage(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);
    }
    /** @}*/ // end of sytrd_2stage

    /** @defgroup hetrd_2stage hetrd_2stage
     * @ingroup SYM
     * @{
     */
    /*! @brief HETRD_2STAGE reduces a complex Hermitian matrix A to real symmetric \n
     tridiagonal form T
 * @details
 * \b Purpose:
    \verbatim
    HETRD_2STAGE reduces a complex Hermitian matrix A to real symmetric
    tridiagonal form T by a unitary similarity transformation:
    Q1**H Q2**H* A * Q2 * Q1 = T.
    \endverbatim

 * @param[in] VECT
          VECT is CHARACTER*1 \n
          = 'N':  No need for the Housholder representation,
                  in particular for the second stage (Band to
                  tridiagonal) and thus LHOUS2 is of size fla_max(1, 4*N); \n
          = 'V':  the Householder representation is needed to
                  either generate Q1 Q2 or to apply Q1 Q2,
                  then LHOUS2 is to be queried and computed.
                  (NOT AVAILABLE IN THIS RELEASE). \n
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
          On exit, if UPLO = 'U', the band superdiagonal
          of A are overwritten by the corresponding elements of the
          internal band-diagonal matrix AB, and the elements above
          the KD superdiagonal, with the array TAU, represent the unitary
          matrix Q1 as a product of elementary reflectors; if UPLO
          = 'L', the diagonal and band subdiagonal of A are over-
          written by the corresponding elements of the internal band-diagonal
          matrix AB, and the elements below the KD subdiagonal, with
          the array TAU, represent the unitary matrix Q1 as a product
          of elementary reflectors. See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] D
          D is REAL array, dimension (N) \n
          The diagonal elements of the tridiagonal matrix T. \n
 * @param[out] E
          E is REAL array, dimension (N-1) \n
          The off-diagonal elements of the tridiagonal matrix T. \n
 * @param[out] TAU
          TAU is COMPLEX array, dimension (N-KD) \n
          The scalar factors of the elementary reflectors of
          the first stage (see Further Details). \n
 * @param[out] HOUS2
          HOUS2 is COMPLEX array, dimension (LHOUS2) \n
          Stores the Householder representation of the stage2
          band to tridiagonal. \n
 * @param[in] LHOUS2
          LHOUS2 is INTEGER \n
          The dimension of the array HOUS2. \n
          If LWORK = -1, or LHOUS2=-1,
          then a query is assumed; the routine
          only calculates the optimal size of the HOUS2 array,   returns
          this value as the first entry of the HOUS2 array, and no error
          message related to LHOUS2 is issued by XERBLA. \n
          If VECT='N', LHOUS2 = fla_max(1, 4*n); \n
          if VECT='V', option not yet available. \n
 * @param[out] WORK
          WORK is COMPLEX array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK = MAX(1, dimension) \n
          If LWORK = -1, or LHOUS2 = -1,
          then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array,   returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
          LWORK = MAX(1, dimension) where \n
          dimension   = fla_max(stage1,stage2) + (KD+1)*N \n
                      = N*KD + N*max(KD+1,FACTOPTNB)  \n
                        + fla_max(2*KD*KD, KD*NTHREADS)  \n
                        + (KD+1)*N  \n
          where KD is the blocking size of the reduction,
          FACTOPTNB is the blocking used by the QR or LQ
          algorithm, usually FACTOPTNB=128 is a good choice
          NTHREADS is the number of threads used when
          openMP compilation is enabled, otherwise =1. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T, typename Ta>
    void hetrd_2stage(char *vect, char *uplo, integer *n, T *a, integer *lda, Ta *d, Ta *e, T *tau,
                      T *hous2, integer *lhous2, T *work, integer *lwork, integer *info)
    {
        hetrd_2stage(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);
    }
    /** @}*/ // end of hetrd_2stage

    /** @defgroup hetrd_he2hb hetrd_he2hb
     * @ingroup SYM
     * @{
     */
    /*! @brief HETRD_HE2HB reduces a complex Hermitian matrix A to complex Hermitian
 band-diagonal form AB

 * @details
 * \b Purpose:
    \verbatim
    HETRD_HE2HB reduces a complex Hermitian matrix A to complex Hermitian
    band-diagonal form AB by a unitary similarity transformation:
    Q**H * A * Q = AB.
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
           The number of superdiagonals of the reduced matrix if UPLO = 'U',
           or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
           The reduced matrix is stored in the array AB. \n
  * @param[in,out] A
           A is COMPLEX array, dimension (LDA,N) \n
           On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
           N-by-N upper triangular part of A contains the upper
           triangular part of the matrix A, and the strictly lower
           triangular part of A is not referenced.  If UPLO = 'L', the
           leading N-by-N lower triangular part of A contains the lower
           triangular part of the matrix A, and the strictly upper
           triangular part of A is not referenced. \n \n
           On exit, if UPLO = 'U', the diagonal and first superdiagonal
           of A are overwritten by the corresponding elements of the
           tridiagonal matrix T, and the elements above the first
           superdiagonal, with the array TAU, represent the unitary
           matrix Q as a product of elementary reflectors; if UPLO
           = 'L', the diagonal and first subdiagonal of A are over-
           written by the corresponding elements of the tridiagonal
           matrix T, and the elements below the first subdiagonal, with
           the array TAU, represent the unitary matrix Q as a product
           of elementary reflectors. See Further Details. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,N). \n
  * @param[out] AB
           AB is COMPLEX array, dimension (LDAB,N) \n
           On exit, the upper or lower triangle of the Hermitian band
           matrix A, stored in the first KD+1 rows of the array.  The
           j-th column of A is stored in the j-th column of the array AB
           as follows: \n
           if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
           if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
  * @param[in] LDAB
           LDAB is INTEGER \n
           The leading dimension of the array AB.  LDAB >= KD+1. \n
  * @param[out] TAU
           TAU is COMPLEX array, dimension (N-KD) \n
           The scalar factors of the elementary reflectors (see Further
           Details). \n
  * @param[out] WORK
           WORK is COMPLEX array, dimension (LWORK) \n
           On exit, if INFO = 0, or if LWORK=-1,
           WORK(1)   returns the size of LWORK. \n
  * @param[in] LWORK
           LWORK is INTEGER \n
           The dimension of the array WORK which should be calculated
           by a workspace query. LWORK = MAX(1, LWORK_QUERY) \n
           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array,   returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA. \n
           LWORK_QUERY = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD \n
           where FACTOPTNB is the blocking used by the QR or LQ
           algorithm, usually FACTOPTNB=128 is a good choice otherwise
           putting LWORK=-1 will provide the size of WORK. \n
  * @param[out] INFO
           INFO is INTEGER \n
           = 0:  successful exit \n
           < 0:  if INFO = -i, the i-th argument had an illegal value  \n

 *  * */
    template <typename T>
    void hetrd_he2hb(char *uplo, integer *n, integer *kd, T *a, integer *lda, T *ab, integer *ldab,
                     T *tau, T *work, integer *lwork, integer *info)
    {
        hetrd_he2hb(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
    }
    /** @}*/ // end of hetrd_he2hb

    /** @defgroup sytrd_hb2st sytrd_hb2st
     * @ingroup SYM
     * @{
     */
    /*! @brief SYTRD_SB2ST reduces a real symmetric band matrix A to real symmetric tridiagonal
 form T

 * @details
 * \b Purpose:
    \verbatim
    SYTRD_SB2ST reduces a real symmetric band matrix A to real symmetric
    tridiagonal form T by a orthogonal similarity transformation:
    Q**T * A * Q = T.
    \endverbatim

  * @param[in] STAGE1
           STAGE1 is CHARACTER*1 \n
           = 'N':  "No": to mention that the stage 1 of the reduction
                   from dense to band using the ssytrd_sy2sb routine
                   was not called before this routine to reproduce AB.
                   In other term this routine is called as standalone.  \n
           = 'Y':  "Yes": to mention that the stage 1 of the
                   reduction from dense to band using the ssytrd_sy2sb
                   routine has been called to produce AB (e.g., AB is
                   the output of ssytrd_sy2sb. \n
  * @param[in] VECT
           VECT is CHARACTER*1 \n
           = 'N':  No need for the Housholder representation,
                   and thus LHOUS is of size fla_max(1, 4*N); \n
           = 'V':  the Householder representation is needed to
                   either generate or to apply Q later on,
                   then LHOUS is to be queried and computed.
                   (NOT AVAILABLE IN THIS RELEASE). \n
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
           On exit, the diagonal elements of AB are overwritten by the
           diagonal elements of the tridiagonal matrix T; if KD > 0, the
           elements on the first superdiagonal (if UPLO = 'U') or the
           first subdiagonal (if UPLO = 'L') are overwritten by the
           off-diagonal elements of T; the rest of AB is overwritten by
           values generated during the reduction. \n
  * @param[in] LDAB
           LDAB is INTEGER \n
           The leading dimension of the array AB.  LDAB >= KD+1. \n
  * @param[out] D
           D is REAL array, dimension (N) \n
           The diagonal elements of the tridiagonal matrix T. \n
  * @param[out] E
           E is REAL array, dimension (N-1) \n
           The off-diagonal elements of the tridiagonal matrix T: \n
           E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'. \n
  * @param[out] HOUS
           HOUS is REAL array, dimension LHOUS, that
           store the Householder representation. \n
  * @param[in] LHOUS
           LHOUS is INTEGER \n
           The dimension of the array HOUS. LHOUS = MAX(1, dimension) \n
           If LWORK = -1, or LHOUS=-1,
           then a query is assumed; the routine
           only calculates the optimal size of the HOUS array,   returns
           this value as the first entry of the HOUS array, and no error
           message related to LHOUS is issued by XERBLA. \n
           LHOUS = MAX(1, dimension) where \n
           dimension = 4*N if VECT='N' \n
           not available now if VECT='H' \n
  * @param[out] WORK
           WORK is REAL array, dimension LWORK. \n
  * @param[in] LWORK
           LWORK is INTEGER \n
           The dimension of the array WORK. LWORK = MAX(1, dimension) \n
           If LWORK = -1, or LHOUS=-1,
           then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array,   returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA. \n
           LWORK = MAX(1, dimension) where
           dimension   = (2KD+1)*N + KD*NTHREADS \n
           where KD is the blocking size of the reduction,
           FACTOPTNB is the blocking used by the QR or LQ
           algorithm, usually FACTOPTNB=128 is a good choice
           NTHREADS is the number of threads used when
           openMP compilation is enabled, otherwise =1. \n
  * @param[out] INFO
           INFO is INTEGER \n
           = 0:  successful exit \n
           < 0:  if INFO = -i, the i-th argument had an illegal value  \n

 *  * */
    template <typename T>
    void sytrd_sb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, T *ab,
                     integer *ldab, T *d, T *e, T *hous, integer *lhous, T *work, integer *lwork,
                     integer *info)
    {
        hetrd_hb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
    }
    /** @}*/ // end of sytrd_hb2st

    /** @defgroup hb2st_kernels hb2st_kernels
     * @ingroup SYM
     * @{
     */
    /*! @brief HB2ST_KERNELS is an internal routine used by the HETRD_HB2ST subroutine

 * @details
 * \b Purpose:
    \verbatim
    HB2ST_KERNELS is an internal routine used by the CHETRD_HB2ST
    subroutine.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
 * @param[in] WANTZ
          WANTZ is LOGICAL which indicate if Eigenvalue are requested or both
          Eigenvalue/Eigenvectors. \n
 * @param[in] TTYPE
          TTYPE is INTEGER \n
 * @param[in] ST
          ST is INTEGER \n
          internal parameter for indices. \n
 * @param[in] ED
          ED is INTEGER \n
          internal parameter for indices. \n
 * @param[in] SWEEP
          SWEEP is INTEGER \n
          internal parameter for indices. \n
 * @param[in] N
          N is INTEGER. The order of the matrix A. \n
 * @param[in] NB
          NB is INTEGER. The size of the band. \n
 * @param[in] IB
          IB is INTEGER. \n
 * @param[in, out] A
          A is COMPLEX array. A pointer to the matrix A. \n
 * @param[in] LDA
          LDA is INTEGER. The leading dimension of the matrix A. \n
 * @param[out] V
          V is COMPLEX array, dimension 2*n if eigenvalues only are
          requested or to be queried for vectors. \n
 * @param[out] TAU
          TAU is COMPLEX array, dimension (2*n).
          The scalar factors of the Householder reflectors are stored
          in this array. \n
 * @param[in] LDVT
          LDVT is INTEGER. \n
 * @param[out] WORK
          WORK is COMPLEX array. Workspace of size nb. \n

 *  * */
    template <typename T>
    void hb2st_kernels(char *uplo, logical *wantz, integer *ttype, integer *st, integer *ed,
                       integer *sweep, integer *n, integer *nb, integer *ib, T *a, integer *lda,
                       T *v, T *tau, integer *ldvt, T *work)
    {
        hb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
    }
    /** @}*/ // end of hb2st_kernels

    /** @defgroup sb2st_kernels sb2st_kernels
     * @ingroup SYM
     * @{
     */
    /*! @brief SB2ST_KERNELS is an internal routine used by the SYTRD_SB2ST subroutine

 * @details
 * \b Purpose:
    \verbatim
    SB2ST_KERNELS is an internal routine used by the SYTRD_SB2ST
    subroutine.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
 * @param[in] WANTZ
          WANTZ is LOGICAL which indicate if Eigenvalue are requested or both
          Eigenvalue/Eigenvectors. \n
 * @param[in] TTYPE
          TTYPE is INTEGER \n
 * @param[in] ST
          ST is INTEGER \n
          internal parameter for indices. \n
 * @param[in] ED
          ED is INTEGER \n
          internal parameter for indices. \n
 * @param[in] SWEEP
          SWEEP is INTEGER \n
          internal parameter for indices. \n
 * @param[in] N
          N is INTEGER. The order of the matrix A. \n
 * @param[in] NB
          NB is INTEGER. The size of the band. \n
 * @param[in] IB
          IB is INTEGER. \n
 * @param[in, out] A
          A is REAL array. A pointer to the matrix A. \n
 * @param[in] LDA
          LDA is INTEGER. The leading dimension of the matrix A. \n
 * @param[out] V
          V is REAL array, dimension 2*n if eigenvalues only are
          requested or to be queried for vectors. \n
 * @param[out] TAU
          TAU is REAL array, dimension (2*n). \n
          The scalar factors of the Householder reflectors are stored
          in this array. \n
 * @param[in] LDVT
          LDVT is INTEGER. \n
 * @param[out] WORK
          WORK is REAL array. Workspace of size nb. \n

 *  * */
    template <typename T>
    void sb2st_kernels(char *uplo, logical *wantz, integer *ttype, integer *st, integer *ed,
                       integer *sweep, integer *n, integer *nb, integer *ib, T *a, integer *lda,
                       T *v, T *tau, integer *ldvt, T *work)
    {
        sb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
    }
    /** @}*/ // end of sb2st_kernel

    /** @defgroup lae2 lae2
     * @ingroup SYM
     * @{
     */
    /*! @brief LAE2 computes the eigenvalues of a 2-by-2 symmetric matrix

 * @details
 * \b Purpose:
    \verbatim
    LAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
       [  A   B  ]
       [  B   C  ].
    On   return, RT1 is the eigenvalue of larger absolute value, and RT2
    is the eigenvalue of smaller absolute value.
    \endverbatim

 * @param[in] A
          A is REAL \n
          The (1,1) element of the 2-by-2 matrix. \n
 * @param[in] B
          B is REAL \n
          The (1,2) and (2,1) elements of the 2-by-2 matrix. \n
 * @param[in] C
          C is REAL \n
          The (2,2) element of the 2-by-2 matrix. \n
 * @param[out] RT1
          RT1 is REAL \n
          The eigenvalue of larger absolute value. \n
 * @param[out] RT2
          RT2 is REAL \n
          The eigenvalue of smaller absolute value. \n

 *  * */
    template <typename T>
    void lae2(T *a, T *b, T *c, T *rt1, T *rt2)
    {
        lae2(a, b, c, rt1, rt2);
    }
    /** @}*/ // end of lae2

    /** @defgroup laesy laesy
     * @ingroup SYM
     * @{
     */
    /*! @brief LAESY computes the eigenvalues and eigenvectors of a 2-by-2 complex symmetric
 matrix

 * @details
 * \b Purpose:
    \verbatim
    LAESY computes the eigendecomposition of a 2-by-2 symmetric matrix
       ((A, B);(B, C))
    provided the norm of the matrix of eigenvectors is larger than
    some threshold value.

    RT1 is the eigenvalue of larger absolute value, and RT2 of
    smaller absolute value.  If the eigenvectors are computed, then
    on   return (CS1, SN1) is the unit eigenvector for RT1, hence

    [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]
    [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]
    \endverbatim

 * @param[in] A
          A is COMPLEX \n
          The (1, 1) element of input matrix. \n
 * @param[in] B
          B is COMPLEX \n
          The (1, 2) element of input matrix.  The (2, 1) element
          is also given by B, since the 2-by-2 matrix is symmetric. \n
 * @param[in] C
          C is COMPLEX \n
          The (2, 2) element of input matrix. \n
 * @param[out] RT1
          RT1 is COMPLEX \n
          The eigenvalue of larger modulus. \n
 * @param[out] RT2
          RT2 is COMPLEX \n
          The eigenvalue of smaller modulus. \n
 * @param[out] EVSCAL
          EVSCAL is COMPLEX \n
          The complex value by which the eigenvector matrix was scaled
          to make it orthonormal.  If EVSCAL is zero, the eigenvectors
          were not computed.  This means one of two things:  the 2-by-2
          matrix could not be diagonalized, or the norm of the matrix
          of eigenvectors before scaling was larger than the threshold
          value THRESH (set below). \n
 * @param[out] CS1
          CS1 is COMPLEX \n
 * @param[out] SN1
          SN1 is COMPLEX \n
          If EVSCAL .NE. 0, (CS1, SN1) is the unit right eigenvector
          for RT1. \n

 * */
    template <typename T>
    void laesy(T *a, T *b, T *c, T *rt1, T *rt2, T *evscal, T *cs1, T *sn1)
    {
        laesy(a, b, c, rt1, rt2, evscal, cs1, sn1);
    }
    /** @}*/ // end of laesy

    /** @defgroup laev2 laev2
     * @ingroup SYM
     * @{
     */
    /*! @brief LAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian
 matrix

 * @details
 * \b Purpose:
    \verbatim
    LAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
       [  A   B  ]
       [  B   C  ].
    On   return, RT1 is the eigenvalue of larger absolute value, RT2 is the
    eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
    eigenvector for RT1, giving the decomposition

       [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
       [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
    \endverbatim

 * @param[in] A
          A is REAL \n
          The (1,1) element of the 2-by-2 matrix. \n
 * @param[in] B
          B is REAL \n
          The (1,2) element and the conjugate of the (2,1) element of
          the 2-by-2 matrix. \n
 * @param[in] C
          C is REAL \n
          The (2,2) element of the 2-by-2 matrix. \n
 * @param[out] RT1
          RT1 is REAL \n
          The eigenvalue of larger absolute value. \n
 * @param[out] RT2
          RT2 is REAL \n
          The eigenvalue of smaller absolute value. \n
 * @param[out] CS1
          CS1 is REAL \n
 * @param[out] SN1
          SN1 is REAL \n
          The vector (CS1, SN1) is a unit right eigenvector for RT1.  \n

 * */
    template <typename T>
    void laev2(T *a, T *b, T *c, T *rt1, T *rt2, T *cs1, T *sn1)
    {
        laev2(a, b, c, rt1, rt2, cs1, sn1);
    }
    template <typename T, typename Ta>
    void laev2(T *a, T *b, T *c, Ta *rt1, Ta *rt2, Ta *cs1, T *sn1)
    {
        laev2(a, b, c, rt1, rt2, cs1, sn1);
    }
    /** @}*/ // end of laev2

    /** @defgroup lagtf lagtf
     * @ingroup SYM
     * @{
     */
    /*! @brief LAGTF computes an LU factorization of a matrix T-Î»I, where T is a general  \n
     tridiagonal matrix, and Î» a scalar, using partial pivoting with row interchanges
 * @details
 * \b Purpose:
    \verbatim
    LAGTF factorizes the matrix (T - lambda*I), where T is an n by n
    tridiagonal matrix and lambda is a scalar, as

       T - lambda*I = PLU,

    where P is a permutation matrix, L is a unit lower tridiagonal matrix
    with at most one non-zero sub-diagonal elements per column and U is
    an upper triangular matrix with at most two non-zero super-diagonal
    elements per column.

    The factorization is obtained by Gaussian elimination with partial
    pivoting and implicit row scaling.

    The parameter LAMBDA is included in the routine so that SLAGTF may
    be used, in conjunction with SLAGTS, to obtain eigenvectors of T by
    inverse iteration.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. \n
 * @param[in,out] A
          A is REAL array, dimension (N) \n
          On entry, A must contain the diagonal elements of T.
 \n
          On exit, A is overwritten by the n diagonal elements of the
          upper triangular matrix U of the factorization of T. \n
 * @param[in] LAMBDA
          LAMBDA is REAL \n
          On entry, the scalar lambda. \n
 * @param[in,out] B
          B is REAL array, dimension (N-1) \n
          On entry, B must contain the (n-1) super-diagonal elements of
          T.
 \n
          On exit, B is overwritten by the (n-1) super-diagonal
          elements of the matrix U of the factorization of T. \n
 * @param[in,out] C
          C is REAL array, dimension (N-1) \n
          On entry, C must contain the (n-1) sub-diagonal elements of
          T.
 \n
          On exit, C is overwritten by the (n-1) sub-diagonal elements
          of the matrix L of the factorization of T. \n
 * @param[in] TOL
          TOL is REAL \n
          On entry, a relative tolerance used to indicate whether or
          not the matrix (T - lambda*I) is nearly singular. TOL should
          normally be chose as approximately the largest relative error
          in the elements of T. For example, if the elements of T are
          correct to about 4 significant figures, then TOL should be
          set to about 5*10**(-4). If TOL is supplied as less than eps,
          where eps is the relative machine precision, then the value
          eps is used in place of TOL. \n
 * @param[out] D
          D is REAL array, dimension (N-2) \n
          On exit, D is overwritten by the (n-2) second super-diagonal
          elements of the matrix U of the factorization of T. \n
 * @param[out] IN
          IN is INTEGER array, dimension (N) \n
          On exit, IN contains details of the permutation matrix P. If
          an interchange occurred at the kth step of the elimination,
          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)
            returns the smallest positive integer j such that
 \n
             abs(u(j,j)) <= norm((T - lambda*I)(j))*TOL,
 \n
          where norm(A(j)) denotes the sum of the absolute values of
          the jth row of the matrix A. If no such j exists then IN(n)
          is   returned as zero. If IN(n) is   returned as positive, then a
          diagonal element of U is small, indicating that
          (T - lambda*I) is singular or nearly singular, \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -k, the kth argument had an illegal value \n

 * */
    template <typename T>
    void lagtf(integer *n, T *a, T *lambda, T *b, T *c, T *tol, T *d, integer *in, integer *info)
    {
        lagtf(n, a, lambda, b, c, tol, d, in, info);
    }
    /** @}*/ // end of lagtf

    /** @defgroup lagts lagts
     * @ingroup SYM
     * @{
     */
    /*! @brief LAGTS solves the system of equations (T-Î»I)x = y or (T-Î»I)Tx = y,where T is a
 general \n tridiagonal matrix and Î» a scalar, using the LU factorization computed by slagtf
 * @details
 * \b Purpose:
    \verbatim
    LAGTS may be used to solve one of the systems of equations

       (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y,

    where T is an n by n tridiagonal matrix, for x, following the
    factorization of (T - lambda*I) as

       (T - lambda*I) = P*L*U ,

    by routine LAGTF. The choice of equation to be solved is
    controlled by the argument JOB, and in each case there is an option
    to perturb zero or very small diagonal elements of U, this option
    being intended for use in applications such as inverse iteration.
    \endverbatim

 * @param[in] JOB
          JOB is INTEGER \n
          Specifies the job to be performed by SLAGTS as follows: \n
          =  1: The equations  (T - lambda*I)x = y  are to be solved,
                but diagonal elements of U are not to be perturbed. \n
          = -1: The equations  (T - lambda*I)x = y  are to be solved
                and, if overflow would otherwise occur, the diagonal
                elements of U are to be perturbed. See argument TOL
                below. \n
          =  2: The equations  (T - lambda*I)**Tx = y  are to be solved,
                but diagonal elements of U are not to be perturbed. \n
          = -2: The equations  (T - lambda*I)**Tx = y  are to be solved
                and, if overflow would otherwise occur, the diagonal
                elements of U are to be perturbed. See argument TOL
                below. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix T. \n
 * @param[in] A
          A is REAL array, dimension (N) \n
          On entry, A must contain the diagonal elements of U as
          returned from SLAGTF. \n
 * @param[in] B
          B is REAL array, dimension (N-1) \n
          On entry, B must contain the first super-diagonal elements of
          U as returned from SLAGTF. \n
 * @param[in] C
          C is REAL array, dimension (N-1) \n
          On entry, C must contain the sub-diagonal elements of L as
          returned from SLAGTF. \n
 * @param[in] D
          D is REAL array, dimension (N-2) \n
          On entry, D must contain the second super-diagonal elements
          of U as returned from SLAGTF. \n
 * @param[in] IN
          IN is INTEGER array, dimension (N) \n
          On entry, IN must contain details of the matrix P as   returned
          from SLAGTF. \n
 * @param[in,out] Y
          Y is REAL array, dimension (N) \n
          On entry, the right hand side vector y.
          On exit, Y is overwritten by the solution vector x. \n
 * @param[in,out] TOL
          TOL is REAL \n
          On entry, with  JOB < 0, TOL should be the minimum
          perturbation to be made to very small diagonal elements of U.
          TOL should normally be chosen as about eps*norm(U), where eps
          is the relative machine precision, but if TOL is supplied as
          non-positive, then it is reset to eps*max(abs(u(i,j))).
          If  JOB > 0  then TOL is not referenced.
 \n
          On exit, TOL is changed as described above, only if TOL is
          non-positive on entry. Otherwise TOL is unchanged. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: overflow would occur when computing the INFO(th)
               element of the solution vector x. This can only occur
               when JOB is supplied as positive and either means
               that a diagonal element of U is very small, or that
               the elements of the right-hand side vector y are very
               large. \n

 * */
    template <typename T>
    void lagts(integer *job, integer *n, T *a, T *b, T *c, T *d, integer *in, T *y, T *tol,
               integer *info)
    {
        lagts(job, n, a, b, c, d, in, y, tol, info);
    }
    /** @}*/ // end of lagts

    /** @defgroup sptrd sptrd
     * @ingroup SYM
     * @{
     */
    /*! @brief SPTRD reduces a real symmetric matrix A stored in packed form to symmetric
 tridiagonal form T

 * @details
 * \b Purpose:
    \verbatim
     SPTRD reduces a real symmetric matrix A stored in packed form to
     symmetric tridiagonal form T by an orthogonal similarity
     transformation: Q**T * A * Q = T.
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
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
          On exit, if UPLO = 'U', the diagonal and first superdiagonal
          of A are overwritten by the corresponding elements of the
          tridiagonal matrix T, and the elements above the first
          superdiagonal, with the array TAU, represent the orthogonal
          matrix Q as a product of elementary reflectors; if UPLO
          = 'L', the diagonal and first subdiagonal of A are over-
          written by the corresponding elements of the tridiagonal
          matrix T, and the elements below the first subdiagonal, with
          the array TAU, represent the orthogonal matrix Q as a product
          of elementary reflectors. See Further Details. \n
 * @param[out] D
          D is REAL array, dimension (N) \n
          The diagonal elements of the tridiagonal matrix T:
          D(i) = A(i,i). \n
 * @param[out] E
          E is REAL array, dimension (N-1) \n
          The off-diagonal elements of the tridiagonal matrix T: \n
          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. \n
 * @param[out] TAU
          TAU is REAL array, dimension (N-1) \n
          The scalar factors of the elementary reflectors (see Further
          Details). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sptrd(char *uplo, integer *n, T *ap, T *d, T *e, T *tau, integer *info)
    {
        sptrd(uplo, n, ap, d, e, tau, info);
    }
    template <typename T, typename Ta>
    void hptrd(char *uplo, integer *n, T *ap, Ta *d, Ta *e, T *tau, integer *info)
    {
        hptrd(uplo, n, ap, d, e, tau, info);
    }
    /** @}*/ // end of sptrd

    /** @defgroup opgtr {up,op}gtr
     * @ingroup SYM
     * @{
     */
    /*! @brief OPGTR generates a real orthogonal matrix Q

 * @details
 * \b Purpose:
    \verbatim
     OPGTR generates a real orthogonal matrix Q which is defined as the
     product of n-1 elementary reflectors H(i) of order n, as returned by
     SSPTRD using packed storage:

     if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),

     if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U': Upper triangular packed storage used in previous
                 call to SSPTRD; \n
          = 'L': Lower triangular packed storage used in previous
                 call to SSPTRD. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix Q. N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The vectors which define the elementary reflectors, as
          returned by SSPTRD. \n
 * @param[in] TAU
          TAU is REAL array, dimension (N-1) \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by SSPTRD. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ,N) \n
          The N-by-N orthogonal matrix Q. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. LDQ >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N-1) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void opgtr(char *uplo, integer *n, T *ap, T *tau, T *q, integer *ldq, T *work, integer *info)
    {
        opgtr(uplo, n, ap, tau, q, ldq, work, info);
    }
    template <typename T>
    void upgtr(char *uplo, integer *n, T *ap, T *tau, T *q, integer *ldq, T *work, integer *info)
    {
        upgtr(uplo, n, ap, tau, q, ldq, work, info);
    }
    /** @}*/ // end of opgtr

    /** @defgroup opmtr {op,up}mtr
     * @ingroup SYM
     * @{
     */
    /*! @brief OPMTR overwrites the general real M-by-N matrix

 * @details
 * \b Purpose:
    \verbatim
     OPMTR overwrites the general real M-by-N matrix C with
                     SIDE = 'L'     SIDE = 'R'
     TRANS = 'N':      Q * C          C * Q
     TRANS = 'T':      Q**T * C       C * Q**T

     where Q is a real orthogonal matrix of order nq, with nq = m if
     SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     nq-1 elementary reflectors, as returned by SSPTRD using packed
     storage:
     if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
     if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
    \endverbatim

 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          = 'L': apply Q or Q**T from the Left; \n
          = 'R': apply Q or Q**T from the Right. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U': Upper triangular packed storage used in previous
                 call to SSPTRD; \n
          = 'L': Lower triangular packed storage used in previous
                 call to SSPTRD. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N':  No transpose, apply Q; \n
          = 'T':  Transpose, apply Q**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix C. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix C. N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension \n
                               (M*(M+1)/2) if SIDE = 'L' \n
                               (N*(N+1)/2) if SIDE = 'R' \n
          The vectors which define the elementary reflectors, as
          returned by SSPTRD.  AP is modified by the routine but
          restored on exit. \n
 * @param[in] TAU
          TAU is REAL array, dimension (M-1) if SIDE = 'L'
                                     or (N-1) if SIDE = 'R' \n
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by SSPTRD. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,N) \n
          On entry, the M-by-N matrix C. \n
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >= fla_max(1,M). \n
 * @param[out]	WORK
          WORK is REAL array, dimension \n
                                   (N) if SIDE = 'L' \n
                                   (M) if SIDE = 'R' \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void opmtr(char *side, char *uplo, char *trans, integer *m, integer *n, T *ap, T *tau, T *c,
               integer *ldc, T *work, integer *info)
    {
        opmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
    }

    template <typename T>
    void upmtr(char *side, char *uplo, char *trans, integer *m, integer *n, T *ap, T *tau, T *c,
               integer *ldc, T *work, integer *info)
    {
        upmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
    }
    /** @}*/ // end of opmtr

    /** @defgroup sbtrd {hb,sb}trd
     * @ingroup SYM
     * @{
     */
    /*! @brief SBTRD reduces a real symmetric band matrix A to symmetric tridiagonal form T

 * @details
 * \b Purpose:
    \verbatim
     SBTRD reduces a real symmetric band matrix A to symmetric
     tridiagonal form T by an orthogonal similarity transformation:
     Q**T * A * Q = T.
    \endverbatim

 * @param[in] VECT
          VECT is CHARACTER*1 \n
          = 'N':  do not form Q; \n
          = 'V':  form Q; \n
          = 'U':  update a matrix X, by forming X*Q. \n
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
          On exit, the diagonal elements of AB are overwritten by the
          diagonal elements of the tridiagonal matrix T; if KD > 0, the
          elements on the first superdiagonal (if UPLO = 'U') or the
          first subdiagonal (if UPLO = 'L') are overwritten by the
          off-diagonal elements of T; the rest of AB is overwritten by
          values generated during the reduction. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[out] D
          D is REAL array, dimension (N) \n
          The diagonal elements of the tridiagonal matrix T. \n
 * @param[out] E
          E is REAL array, dimension (N-1) \n
          The off-diagonal elements of the tridiagonal matrix T: \n
          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if VECT = 'U', then Q must contain an N-by-N
          matrix X; if VECT = 'N' or 'V', then Q need not be set. \n
 \n
          On exit: \n
          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q; \n
          if VECT = 'U', Q contains the product X*Q; \n
          if VECT = 'N', the array Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. \n
          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sbtrd(char *vect, char *uplo, integer *n, integer *kd, T *ab, integer *ldab, T *d, T *e,
               T *q, integer *ldq, T *work, integer *info)
    {
        sbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    }
    template <typename T, typename Ta>
    void hbtrd(char *vect, char *uplo, integer *n, integer *kd, T *ab, integer *ldab, Ta *d, Ta *e,
               T *q, integer *ldq, T *work, integer *info)
    {
        hbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    }
    /** @}*/ // end of sbtrd

    /** @defgroup sygst sygst
     * @ingroup SYM
     * @{
     */
    /*! @brief Reduction of a real symmetric-definite generalized eigenproblem to standard form
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric-definite generalized eigenproblem to standard form.

        If itype = 1, the problem is A*X = lambda*B*X,
        and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**T or L**T*A*L.

        B must have been previously factorized as U**T*U or L*L**T by SPOTRF.
    \endverbatim

    * @param[in] itype
              itype is integer* \n
              = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); \n
              = 2 or 3: compute U*A*U**T or L**T*A*L. \n
    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored and b is factored as
              U**T*U; \n
              = 'L':  Lower triangle of a is stored and b is factored as
              L*L**T. \n
    * @param[in] n
              n is integer* \n
              The order of the matrices a and b.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the leading
              n-by-n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced.  If uplo = 'L', the
              leading n-by-n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if info = 0, the transformed matrix, stored in the
              same format as a. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[in] b
              b is float/double array, dimension (ldb,n) \n
              The triangular factor from the Cholesky factorization of b,
              as returned by SPOTRF. \n
    * @param[in] ldb
              ldb is integer* \n
              The leading dimension of the array b.  ldb >= fla_max(1,n). \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void sygst(integer *itype, char *uplo, integer *n, T *a, integer *lda, T *b, integer *ldb,
               integer *info)
    {
        sygst(itype, uplo, n, a, lda, b, ldb, info);
    }
    /** @}*/ // end of sygst

    /** @defgroup hegst hegst
     * @ingroup SYM
     * @{
     */ /*! @brief Reduction of a complex Hermitian-definite generalized eigenproblem to standard form
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric-definite generalized eigenproblem to standard form.

        If itype = 1, the problem is A*X = lambda*B*X,
        and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**H or L**H*A*L.

        B must have been previously factorized as U**H*U or L*L**H by SPOTRF.
    \endverbatim

    * @param[in] itype
              itype is integer* \n
              = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H); \n
              = 2 or 3: compute U*A*U**H or L**H*A*L. \n
    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored and b is factored as
              U**H*U; \n
              = 'L':  Lower triangle of a is stored and b is factored as
              L*L**H. \n
    * @param[in] n
              n is integer* \n
              The order of the matrices a and b.  n >= 0. \n
    * @param[in,out] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the leading
              n-by-n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced.  If uplo = 'L', the
              leading n-by-n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if info = 0, the transformed matrix, stored in the
              same format as a.
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[in] b
              b is COMPLEX/COMPLEX*16 array, dimension (ldb,n) \n
              The triangular factor from the Cholesky factorization of b,
              as returned by SPOTRF. \n
    * @param[in] ldb
              ldb is integer* \n
              The leading dimension of the array b.  ldb >= fla_max(1,n). \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void hegst(integer *itype, char *uplo, integer *n, T *a, integer *lda, T *b, integer *ldb,
               integer *info)
    {
        hegst(itype, uplo, n, a, lda, b, ldb, info);
    }
    /** @}*/ // end of hegst

    /** @defgroup hegs2 hegs2
     * @ingroup SYM
     * @{
     */
    /*! @brief Reduction of a Hermitian-definite generalized eigenproblem to standard form
    (unblocked algorithm)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a Hermitian definite generalized eigenproblem to standard form, using the
        factorization results obtained from cpotrf (unblocked algorithm).
        If itype = 1, the problem is A*X = lambda*B*X,
        and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**H or L**H *A*L.

        B must have been previously factorized as U**H *U or L*L**H by CPOTRF.
    \endverbatim

    * @param[in] itype
              itype is integer* \n
              = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H); \n
              = 2 or 3: compute U*A*U**H or L**H *A*L. \n
    * @param[in] uplo
              uplo is char* \n
              Specifies whether the upper or lower triangular part of the
              Hermitian matrix a is stored, and how b has been factorized. \n
              = 'U':  Upper triangular \n
              = 'L':  Lower triangular \n
    * @param[in] n
              n is integer* \n
              The order of the matrices a and b.  n >= 0. \n
    * @param[in,out] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the leading
              n by n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced.  If uplo = 'L', the
              leading n by n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if info = 0, the transformed matrix, stored in the
              same format as a. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[in] b
              b is COMPLEX/COMPLEX*16 array, dimension (ldb,n) \n
              The triangular factor from the Cholesky factorization of b,
              as returned by CPOTRF. \n
    * @param[in] ldb
              ldb is integer* \n
              The leading dimension of the array b.  ldb >= fla_max(1,n). \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit. \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void hegs2(integer *itype, char *uplo, integer *n, T *a, integer *lda, T *b, integer *ldb,
               integer *info)
    {
        hegs2(itype, uplo, n, a, lda, b, ldb, info);
    }
    /** @}*/ // end of hegs2

    /** @defgroup sygs2 sygs2
     * @ingroup SYM
     * @{
     */
    /*! @brief Reduction of a symmetric-definite generalized eigenproblem to standard form
    (unblocked algorithm)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric definite generalized eigenproblem to standard form, using
    the factorization results obtained from spotrf (unblocked algorithm). If itype = 1, the
    problem is A*X = lambda*B*X, and A is overwritten by inv(U**T)*A*inv(U) or
    inv(L)*A*inv(L**T)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**T or L**T *A*L.

        B must have been previously factorized as U**T *U or L*L**T by SPOTRF.
    \endverbatim

    * @param[in] itype
              itype is integer* \n
              = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); \n
              = 2 or 3: compute U*A*U**T or L**T *A*L. \n
    * @param[in] uplo
              uplo is char* \n
              Specifies whether the upper or lower triangular part of the
              symmetric matrix a is stored, and how b has been factorized. \n
              = 'U':  Upper triangular \n
              = 'L':  Lower triangular \n
    * @param[in] n
              n is integer* \n
              The order of the matrices a and b.  n >= 0. \n
    * @param[in,out] a
              a is float/doublearray, dimension (lda,n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the leading
              n by n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced.  If uplo = 'L', the
              leading n by n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if info = 0, the transformed matrix, stored in the
              same format as a. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,n). \n
    * @param[in] b
              b is float/double array, dimension (ldb,n) \n
              The triangular factor from the Cholesky factorization of b,
              as returned by SPOTRF. \n
    * @param[in] ldb
              ldb is integer* \n
              The leading dimension of the array b.  ldb >= fla_max(1,n). \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit. \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void sygs2(integer *itype, char *uplo, integer *n, T *a, integer *lda, T *b, integer *ldb,
               integer *info)
    {
        sygs2(itype, uplo, n, a, lda, b, ldb, info);
    }
    /** @}*/ // end of sygs2

    /** @defgroup spgst spgst
     * @ingroup SYM
     * @{
     */
    /*! @brief SPGST reduces a real symmetric-definite generalized eigenproblem to standard form

 * @details
 * \b Purpose:
    \verbatim
     SPGST reduces a real symmetric-definite generalized eigenproblem
     to standard form, using packed storage.

     If ITYPE = 1, the problem is A*x = lambda*B*x,
     and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)

     If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
     B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.

     B must have been previously factorized as U**T*U or L*L**T by SPPTRF.
    \endverbatim

 * @param[in] ITYPE
          ITYPE is INTEGER \n
          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); \n
          = 2 or 3: compute U*A*U**T or L**T*A*L. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored and B is factored as
                  U**T*U; \n
          = 'L':  Lower triangle of A is stored and B is factored as
                  L*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices A and B.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, if INFO = 0, the transformed matrix, stored in the
          same format as A. \n
 * @param[in] BP
          BP is REAL array, dimension (N*(N+1)/2) \n
          The triangular factor from the Cholesky factorization of B,
          stored in the same format as A, as returned by SPPTRF. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void spgst(integer *itype, char *uplo, integer *n, T *ap, T *bp, integer *info)
    {
        spgst(itype, uplo, n, ap, bp, info);
    }
    template <typename T>
    void hpgst(integer *itype, char *uplo, integer *n, T *ap, T *bp, integer *info)
    {
        hpgst(itype, uplo, n, ap, bp, info);
    }
    /** @}*/ // end of spgst

    /** @defgroup pbstf pbstf
     * @ingroup SYM
     * @{
     */
    /*! @brief PBSTF computes a split Cholesky factorization of a real symmetric \n
     positive definite band matrix

 * @details
 * \b Purpose:
    \verbatim
     PBSTF computes a split Cholesky factorization of a real
     symmetric positive definite band matrix A.

     This routine is designed to be used in conjunction with SSBGST.

     The factorization has the form  A = S**T*S  where S is a band matrix
     of the same bandwidth as A and the following structure:

       S = ( U   )
           ( M  L)

     where U is upper triangular of order m = (n+kd)/2, and L is lower
     triangular of order n-m.
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
          matrix A, stored in the first kd+1 rows of the array. The
          j-th column of A is stored in the j-th column of the array AB
          as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
 \n
          On exit, if INFO = 0, the factor S from the split Cholesky
          factorization A = S**T*S. See Further Details. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, the factorization could not be completed,
               because the updated element a(i,i) was negative; the
               matrix A is not positive definite. \n

 *  * */
    template <typename T>
    void pbstf(char *uplo, integer *n, integer *kb, T *bb, integer *ldbb, integer *info)
    {
        pbstf(uplo, n, kb, bb, ldbb, info);
    }
    /** @}*/ // end of pbstf

    /** @defgroup lag2 lag2
     * @ingroup SYM
     * @{
     */
    /*! @brief LAG2 computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, \n
     with scaling as necessary to avoid over-/underflow
 * @details
 * \b Purpose:
    \verbatim
    LAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue
    problem  A - w B, with scaling as necessary to avoid over-/underflow.

    The scaling factor "s" results in a modified eigenvalue equation

        s A - w B

    where  s  is a non-negative scaling factor chosen so that  w, w B,
    and  s A  do not overflow and, if possible, do not underflow, either.
    \endverbatim

 * @param[in] A
          A is REAL array, dimension (LDA, 2) \n
          On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm
          is less than 1/SAFMIN.  Entries less than
          sqrt(SAFMIN)*norm(A) are subject to being treated as zero. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= 2. \n
 * @param[in] B
          B is REAL array, dimension (LDB, 2) \n
          On entry, the 2 x 2 upper triangular matrix B.  It is
          assumed that the one-norm of B is less than 1/SAFMIN.  The
          diagonals should be at least sqrt(SAFMIN) times the largest
          element of B (in absolute value); if a diagonal is smaller
          than that, then  +/- sqrt(SAFMIN) will be used instead of
          that diagonal. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= 2. \n
 * @param[in] SAFMIN
          SAFMIN is REAL \n
          The smallest positive number s.t. 1/SAFMIN does not
          overflow.  (This should always be SLAMCH('S') -- it is an
          argument in order to avoid having to call SLAMCH frequently.) \n
 * @param[out] SCALE1
          SCALE1 is REAL \n
          A scaling factor used to avoid over-/underflow in the
          eigenvalue equation which defines the first eigenvalue.  If
          the eigenvalues are complex, then the eigenvalues are
          (WR1  +/-  WI i) / SCALE1  (which may lie outside the
          exponent range of the machine), SCALE1=SCALE2, and SCALE1
          will always be positive.  If the eigenvalues are real, then
          the first (real) eigenvalue is  WR1 / SCALE1 , but this may
          overflow or underflow, and in fact, SCALE1 may be zero or
          less than the underflow threshold if the exact eigenvalue
          is sufficiently large. \n
 * @param[out] SCALE2
          SCALE2 is REAL \n
          A scaling factor used to avoid over-/underflow in the
          eigenvalue equation which defines the second eigenvalue.  If
          the eigenvalues are complex, then SCALE2=SCALE1.  If the
          eigenvalues are real, then the second (real) eigenvalue is
          WR2 / SCALE2 , but this may overflow or underflow, and in
          fact, SCALE2 may be zero or less than the underflow
          threshold if the exact eigenvalue is sufficiently large. \n
 * @param[out] WR1
          WR1 is REAL \n
          If the eigenvalue is real, then WR1 is SCALE1 times the
          eigenvalue closest to the (2,2) element of A B**(-1).  If the
          eigenvalue is complex, then WR1=WR2 is SCALE1 times the real
          part of the eigenvalues. \n
 * @param[out] WR2
          WR2 is REAL \n
          If the eigenvalue is real, then WR2 is SCALE2 times the
          other eigenvalue.  If the eigenvalue is complex, then
          WR1=WR2 is SCALE1 times the real part of the eigenvalues. \n
 * @param[out] WI
          WI is REAL \n
          If the eigenvalue is real, then WI is zero. If the
          eigenvalue is complex, then WI is SCALE1 times the imaginary
          part of the eigenvalues. WI will always be non-negative. \n

 * */
    template <typename T>
    void lag2(T *a, integer *lda, T *b, integer *ldb, T *safmin, T *scale1, T *scale2, T *wr1,
              T *wr2, T *wi)
    {
        lag2(a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi);
    }
    /** @}*/ // end of lag2

    /** @defgroup laebz laebz
     * @ingroup SYM
     * @{
     */
    /*! @brief LAEBZ computes the number of eigenvalues of a real symmetric tridiagonal \n
     matrix which are less than or equal to a given value, and performs other tasks \n
     required by the routine stebz
 * @details
 * \b Purpose:
    \verbatim
    LAEBZ contains the iteration loops which compute and use the
    function N(w), which is the count of eigenvalues of a symmetric
    tridiagonal matrix T less than or equal to its argument  w.  It
    performs a choice of two types of loops:

    IJOB=1, followed by
    IJOB=2: It takes as input a list of intervals and   returns a list of
            sufficiently small intervals whose union contains the same
            eigenvalues as the union of the original intervals.
            The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
            The output interval (AB(j,1),AB(j,2)] will contain
            eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.

    IJOB=3: It performs a binary search in each input interval
            (AB(j,1),AB(j,2)] for a point  w(j)  such that
            N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
            the search.  If such a w(j) is found, then on output
            AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
            (AB(j,1),AB(j,2)] will be a small interval containing the
            point where N(w) jumps through NVAL(j), unless that point
            lies outside the initial interval.

    Note that the intervals are in all cases half-open intervals,
    i.e., of the form  (a,b] , which includes  b  but not  a .

    To avoid underflow, the matrix should be scaled so that its largest
    element is no greater than  overflow**(1/2) * underflow**(1/4)
    in absolute value.  To assure the most accurate computation
    of small eigenvalues, the matrix should be scaled to be
    not much smaller than that, either.

    See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
    Matrix", Report CS41, Computer Science Dept., Stanford
    University, July 21, 1966

    Note: the arguments are, in general, *not* checked for unreasonable
    values.
    \endverbatim

 * @param[in] IJOB
          IJOB is INTEGER \n
          Specifies what is to be done: \n
          = 1:  Compute NAB for the initial intervals. \n
          = 2:  Perform bisection iteration to find eigenvalues of T. \n
          = 3:  Perform bisection iteration to invert N(w), i.e.,
                to find a point which has a specified number of
                eigenvalues of T to its left. \n
          Other values will cause SLAEBZ to   return with INFO=-1. \n
 * @param[in] NITMAX
          NITMAX is INTEGER \n
          The maximum number of "levels" of bisection to be
          performed, i.e., an interval of width W will not be made
          smaller than 2^(-NITMAX) * W.  If not all intervals
          have converged after NITMAX iterations, then INFO is set
          to the number of non-converged intervals. \n
 * @param[in] N
          N is INTEGER \n
          The dimension n of the tridiagonal matrix T.  It must be at
          least 1. \n
 * @param[in] MMAX
          MMAX is INTEGER \n
          The maximum number of intervals. If more than MMAX intervals
          are generated, then SLAEBZ will quit with INFO=MMAX+1. \n
 * @param[in] MINP
          MINP is INTEGER \n
          The initial number of intervals.  It may not be greater than
          MMAX. \n
 * @param[in] NBMIN
          NBMIN is INTEGER \n
          The smallest number of intervals that should be processed
          using a vector loop.  If zero, then only the scalar loop
          will be used. \n
 * @param[in] ABSTOL
          ABSTOL is REAL \n
          The minimum (absolute) width of an interval.  When an
          interval is narrower than ABSTOL, or than RELTOL times the
          larger (in magnitude) endpoint, then it is considered to be
          sufficiently small, i.e., converged.  This must be at least
          zero. \n
 * @param[in] RELTOL
          RELTOL is REAL \n
          The minimum relative width of an interval.  When an interval
          is narrower than ABSTOL, or than RELTOL times the larger (in
          magnitude) endpoint, then it is considered to be
          sufficiently small, i.e., converged.  Note: this should
          always be at least radix*machine epsilon. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum absolute value of a "pivot" in the Sturm
          sequence loop.
          This must be at least  max |e(j)**2|*safe_min  and at
          least safe_min, where safe_min is at least
          the smallest number that can divide one without overflow. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The diagonal elements of the tridiagonal matrix T. \n
 * @param[in] E
          E is REAL array, dimension (N) \n
          The offdiagonal elements of the tridiagonal matrix T in
          positions 1 through N-1.  E(N) is arbitrary. \n
 * @param[in] E2
          E2 is REAL array, dimension (N) \n
          The squares of the offdiagonal elements of the tridiagonal
          matrix T.  E2(N) is ignored. \n
 * @param[in,out] NVAL
          NVAL is INTEGER array, dimension (MINP) \n
          If IJOB=1 or 2, not referenced. \n
          If IJOB=3, the desired values of N(w).  The elements of NVAL
          will be reordered to correspond with the intervals in AB.
          Thus, NVAL(j) on output will not, in general be the same as
          NVAL(j) on input, but it will correspond with the interval
          (AB(j,1),AB(j,2)] on output. \n
 * @param[in,out] AB
          AB is REAL array, dimension (MMAX,2) \n
          The endpoints of the intervals.  AB(j,1) is  a(j), the left
          endpoint of the j-th interval, and AB(j,2) is b(j), the
          right endpoint of the j-th interval.  The input intervals
          will, in general, be modified, split, and reordered by the
          calculation. \n
 * @param[in,out] C
          C is REAL array, dimension (MMAX) \n
          If IJOB=1, ignored. \n
          If IJOB=2, workspace. \n
          If IJOB=3, then on input C(j) should be initialized to the
          first search point in the binary search. \n
 * @param[out] MOUT
          MOUT is INTEGER \n
          If IJOB=1, the number of eigenvalues in the intervals. \n
          If IJOB=2 or 3, the number of intervals output. \n
          If IJOB=3, MOUT will equal MINP. \n
 * @param[in,out] NAB
          NAB is INTEGER array, dimension (MMAX,2) \n
          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)). \n
          If IJOB=2, then on input, NAB(i,j) should be set.  It must
             satisfy the condition: \n
             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)), \n
             which means that in interval i only eigenvalues
             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,
             NAB(i,j)=N(AB(i,j)), from a previous call to SLAEBZ with
             IJOB=1. \n
             On output, NAB(i,j) will contain
             fla_max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of
             the input interval that the output interval
             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the
             the input values of NAB(k,1) and NAB(k,2). \n
          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),
             unless N(w) > NVAL(i) for all search points  w , in which
             case NAB(i,1) will not be modified, i.e., the output
             value will be the same as the input value (modulo
             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i)
             for all search points  w , in which case NAB(i,2) will
             not be modified.  Normally, NAB should be set to some
             distinctive value(s) before SLAEBZ is called. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MMAX) \n
          Workspace. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (MMAX) \n
          Workspace. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:       All intervals converged. \n
          = 1--MMAX: The last INFO intervals did not converge. \n
          = MMAX+1:  More than MMAX intervals were generated. \n

 * */
    template <typename T>
    void laebz(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp,
               integer *nbmin, T *abstol, T *reltol, T *pivmin, T *d, T *e, T *e2, integer *nval,
               T *ab, T *c, integer *mout, integer *nab, T *work, integer *iwork, integer *info)
    {
        laebz(ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c,
              mout, nab, work, iwork, info);
    }

    /** @}*/ // end of laebz

    /** @defgroup laneg laneg
     * @ingroup SYM
     * @{
     */
    /*! @brief LANEG computes the Sturm count

 * @details
 * \b Purpose:
    \verbatim
    LANEG computes the Sturm count, the number of negative pivots
    encountered while factoring tridiagonal T - sigma I = L D L^T.
    This implementation works directly on the factors without forming
    the tridiagonal matrix T.  The Sturm count is also the number of
    eigenvalues of T less than sigma.

    This routine is called from SLARRB.

    The current routine does not use the PIVMIN parameter but rather
    requires IEEE-754 propagation of Infinities and NaNs.  This
    routine also has no input range restrictions but does require
    default exception handling such that x/0 produces Inf when x is
    non-zero, and Inf/Inf produces NaN.  For more information, see:

      Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in
      Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on
      Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624
      (Tech report version in LAWN 172 with the same title.)
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The N diagonal elements of the diagonal matrix D. \n
 * @param[in] LLD
          LLD is REAL array, dimension (N-1) \n
          The (N-1) elements L(i)*L(i)*D(i). \n
 * @param[in] SIGMA
          SIGMA is REAL \n
          Shift amount in T - sigma I = L D L^T. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot in the Sturm sequence.  May be used
          when zero pivots are encountered on non-IEEE-754
          architectures. \n
 * @param[in] R
          R is INTEGER \n
          The twist index for the twisted factorization that is used
          for the negcount. \n
 * @return INTEGER Return value of the function.
 *  * */
    template <typename T>
    integer laneg(integer *n, T *d, T *lld, T *sigma, T *pivmin, integer *r__)
    {
        return laneg(n, d, lld, sigma, pivmin, r__);
    }
    /** @}*/ // end of laneg

    /** @defgroup laed0 laed0
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED0 used by sstedc. Computes all eigenvalues and corresponding \n
     eigenvectors of an unreduced symmetric tridiagonal matrix using the    \n
     divide and conquer method.
 * @details
 * \b Purpose:
    \verbatim
    LAED0 computes all eigenvalues and corresponding eigenvectors of a
    symmetric tridiagonal matrix using the divide and conquer method.
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          = 0:  Compute eigenvalues only. \n
          = 1:  Compute eigenvectors of original dense symmetric matrix
                also.  On entry, Q contains the orthogonal matrix used
                to reduce the original matrix to tridiagonal form. \n
          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
                matrix. \n
 * @param[in] QSIZ
          QSIZ is INTEGER \n
          The dimension of the orthogonal matrix used to reduce
          the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. \n
 * @param[in] N
          N is INTEGER \n
          The dimension of the symmetric tridiagonal matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the main diagonal of the tridiagonal matrix.
          On exit, its eigenvalues. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The off-diagonal elements of the tridiagonal matrix.
          On exit, E has been destroyed. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ, N) \n
          On entry, Q must contain an N-by-N orthogonal matrix. \n
          If ICOMPQ = 0    Q is not referenced. \n
          If ICOMPQ = 1    On entry, Q is a subset of the columns of the
                           orthogonal matrix used to reduce the full
                           matrix to tridiagonal form corresponding to
                           the subset of the full matrix which is being
                           decomposed at this time. \n
          If ICOMPQ = 2    On entry, Q will be the identity matrix.
                           On exit, Q contains the eigenvectors of the
                           tridiagonal matrix. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  If eigenvectors are
          desired, then  LDQ >= fla_max(1,N).  In any case, LDQ >= 1. \n
 * @param[out] QSTORE
          QSTORE is REAL array, dimension (LDQS, N) \n
          Referenced only when ICOMPQ = 1.  Used to store parts of
          the eigenvector matrix when the updating matrix multiplies
          take place. \n
 * @param[in] LDQS
          LDQS is INTEGER \n
          The leading dimension of the array QSTORE.  If ICOMPQ = 1,
          then  LDQS >= fla_max(1,N).  In any case, LDQS >= 1. \n
 * @param[out] WORK
          WORK is REAL array, \n
          If ICOMPQ = 0 or 1, the dimension of WORK must be at least
                      1 + 3*N + 2*N*lg N + 3*N**2 \n
                      (lg(N) = smallest integer k
                                 such that 2^k >= N) \n
          If ICOMPQ = 2, the dimension of WORK must be at least
                      4*N + N**2. \n
 * @param[out] IWORK
          IWORK is INTEGER array, \n
          If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
                         6 + 6*N + 5*N*lg N. \n
                         (lg(N) = smallest integer k
                                    such that 2^k >= N) \n
          If ICOMPQ = 2, the dimension of IWORK must be at least
                         3 + 5*N. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  The algorithm failed to compute an eigenvalue while
                working on the submatrix lying in rows and columns
                INFO/(N+1) through mod(INFO,N+1).  \n

 * */
    template <typename T>
    void laed0(integer *icompq, integer *qsiz, integer *n, T *d, T *e, T *q, integer *ldq,
               T *qstore, integer *ldqs, T *work, integer *iwork, integer *info)
    {
        laed0(qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info);
    }
    template <typename T, typename Ta>
    void laed0(integer *qsiz, integer *n, Ta *d, Ta *e, T *q, integer *ldq, T *qstore,
               integer *ldqs, Ta *rwork, integer *iwork, integer *info)
    {
        laed0(qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info);
    }
    /** @}*/ // end of laed0

    /** @defgroup laed1 laed1
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED1 used by sstedc. Computes the updated eigensystem of a diagonal \n
     matrix after modification by a rank-one symmetric matrix. Used when the    \n
     original matrix is tridiagonal
 * @details
 * \b Purpose:
    \verbatim
    LAED1 computes the updated eigensystem of a diagonal
    matrix after modification by a rank-one symmetric matrix.  This
    routine is used only for the eigenproblem which requires all
    eigenvalues and eigenvectors of a tridiagonal matrix.  SLAED7 handles
    the case in which eigenvalues only or eigenvalues and eigenvectors
    of a full symmetric matrix (which was reduced to tridiagonal form)
    are desired.

      T = Q(in) (D(in) + RHO * Z*Z**T) Q**T(in) = Q(out) * D(out) * Q**T(out)

       where Z = Q**T*u, u is a vector of length N with ones in the
       CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.

       The eigenvectors of the original matrix are stored in Q, and the
       eigenvalues are in D.  The algorithm consists of three stages:

          The first stage consists of deflating the size of the problem
          when there are multiple eigenvalues or if there is a zero in
          the Z vector.  For each such occurrence the dimension of the
          secular equation problem is reduced by one.  This stage is
          performed by the routine SLAED2.

          The second stage consists of calculating the updated
          eigenvalues. This is done by finding the roots of the secular
          equation via the routine SLAED4 (as called by SLAED3).
          This routine also calculates the eigenvectors of the current
          problem.

          The final stage consists of computing the updated eigenvectors
          directly using the updated eigenvalues.  The eigenvectors for
          the current problem are multiplied with the eigenvectors from
          the overall problem.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The dimension of the symmetric tridiagonal matrix.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the eigenvalues of the rank-1-perturbed matrix.
          On exit, the eigenvalues of the repaired matrix. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, the eigenvectors of the rank-1-perturbed matrix.
          On exit, the eigenvectors of the repaired tridiagonal matrix. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= fla_max(1,N). \n
 * @param[in,out] INDXQ
          INDXQ is INTEGER array, dimension (N) \n
          On entry, the permutation which separately sorts the two
          subproblems in D into ascending order. \n
          On exit, the permutation which will reintegrate the
          subproblems back into sorted order,
          i.e. D(INDXQ(I = 1, N)) will be in ascending order. \n
 * @param[in] RHO
          RHO is REAL \n
          The subdiagonal entry used to create the rank-1 modification. \n
 * @param[in] CUTPNT
          CUTPNT is INTEGER \n
          The location of the last eigenvalue in the leading sub-matrix.
          min(1,N) <= CUTPNT <= N/2. \n
 * @param[out] WORK
          WORK is REAL array, dimension (4*N + N**2) \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (4*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, an eigenvalue did not converge \n

 * */
    template <typename T>
    void laed1(integer *n, T *d, T *q, integer *ldq, integer *indxq, T *rho, integer *cutpnt,
               T *work, integer *iwork, integer *info)
    {
        laed1(*n, d, q, *ldq, *indxq, rho, *cutpnt, work, *iwork, *info);
    }
    /** @}*/ // end of laed1

    /** @defgroup laed2 laed2
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED2 used by sstedc. Merges eigenvalues and deflates secular equation. \n
     Used when the original matrix is tridiagonal
 * @details
 * \b Purpose:
    \verbatim
    LAED2 merges the two sets of eigenvalues together into a single
    sorted set.  Then it tries to deflate the size of the problem.
    There are two ways in which deflation can occur:  when two or more
    eigenvalues are close together or if there is a tiny entry in the
    Z vector.  For each such occurrence the order of the related secular
    equation problem is reduced by one.
    \endverbatim

 * @param[out] K
          K is INTEGER \n
          The number of non-deflated eigenvalues, and the order of the
          related secular equation. 0 <= K <=N. \n
 * @param[in] N
          N is INTEGER \n
          The dimension of the symmetric tridiagonal matrix.  N >= 0. \n
 * @param[in] N1
          N1 is INTEGER \n
          The location of the last eigenvalue in the leading sub-matrix.
          min(1,N) <= N1 <= N/2. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, D contains the eigenvalues of the two submatrices to
          be combined. \n
          On exit, D contains the trailing (N-K) updated eigenvalues
          (those which were deflated) sorted into increasing order. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ, N) \n
          On entry, Q contains the eigenvectors of two submatrices in
          the two square blocks with corners at (1,1), (N1,N1)
          and (N1+1, N1+1), (N,N). \n
          On exit, Q contains the trailing (N-K) updated eigenvectors
          (those which were deflated) in its last N-K columns. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= fla_max(1,N). \n
 * @param[in,out] INDXQ
          INDXQ is INTEGER array, dimension (N) \n
          The permutation which separately sorts the two sub-problems
          in D into ascending order.  Note that elements in the second
          half of this permutation must first have N1 added to their
          values. Destroyed on exit. \n
 * @param[in,out] RHO
          RHO is REAL \n
          On entry, the off-diagonal element associated with the rank-1
          cut which originally split the two submatrices which are now
          being recombined. \n
          On exit, RHO has been modified to the value required by
          SLAED3. \n
 * @param[in] Z
          Z is REAL array, dimension (N) \n
          On entry, Z contains the updating vector (the last
          row of the first sub-eigenvector matrix and the first row of
          the second sub-eigenvector matrix). \n
          On exit, the contents of Z have been destroyed by the updating
          process. \n
 * @param[out] DLAMDA
          DLAMDA is REAL array, dimension (N) \n
          A copy of the first K eigenvalues which will be used by
          SLAED3 to form the secular equation. \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first k values of the final deflation-altered z-vector
          which will be passed to SLAED3. \n
 * @param[out] Q2
          Q2 is REAL array, dimension (N1**2+(N-N1)**2) \n
          A copy of the first K eigenvectors which will be used by
          SLAED3 in a matrix multiply (SGEMM) to solve for the new
          eigenvectors. \n
 * @param[out] INDX
          INDX is INTEGER array, dimension (N) \n
          The permutation used to sort the contents of DLAMDA into
          ascending order. \n
 * @param[out] INDXC
          INDXC is INTEGER array, dimension (N) \n
          The permutation used to arrange the columns of the deflated
          Q matrix into three groups:  the first group contains non-zero
          elements only at and above N1, the second contains
          non-zero elements only below N1, and the third is dense. \n
 * @param[out] INDXP
          INDXP is INTEGER array, dimension (N) \n
          The permutation used to place deflated values of D at the end
          of the array.  INDXP(1:K) points to the nondeflated D-values
          and INDXP(K+1:N) points to the deflated eigenvalues. \n
 * @param[out] COLTYP
          COLTYP is INTEGER array, dimension (N) \n
          During execution, a label which will indicate which of the
          following types a column in the Q2 matrix is: \n
          1 : non-zero in the upper half only; \n
          2 : dense; \n
          3 : non-zero in the lower half only; \n
          4 : deflated. \n
          On exit, COLTYP(i) is the number of columns of type i,
          for i=1 to 4 only. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value.  \n

 * */
    template <typename T>
    void laed2(integer *k, integer *n, integer *n1, T *d, T *q, integer *ldq, integer *indxq,
               T *rho, T *z, T *dlamda, T *w, T *q2, integer *indx, integer *indxc, integer *indxp,
               integer *coltyp, integer *info)
    {
        laed2(*k, *n, *n1, d, q, *ldq, *indxq, rho, z, dlamda, w, q2, *indx, *indxc, *indxp,
              *coltyp, *info);
    }
    /** @}*/ // end of laed2

    /** @defgroup laed3 laed3
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED3 used by sstedc. Finds the roots of the secular equation and \n
     updates the eigenvectors. Used when the original matrix is tridiagonal
 * @details
 * \b Purpose:
    \verbatim
    LAED3 finds the roots of the secular equation, as defined by the
    values in D, W, and RHO, between 1 and K.  It makes the
    appropriate calls to SLAED4 and then updates the eigenvectors by
    multiplying the matrix of eigenvectors of the pair of eigensystems
    being combined by the matrix of eigenvectors of the K-by-K system
    which is solved here.

    This code makes very mild assumptions about floating point
    arithmetic. It will work on machines with a guard digit in
    add/subtract, or on those binary machines without guard digits
    which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
    It could conceivably fail on hexadecimal or decimal machines
    without guard digits, but we know of none.
    \endverbatim

 * @param[in] K
          K is INTEGER \n
          The number of terms in the rational function to be solved by
          SLAED4.  K >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of rows and columns in the Q matrix.
          N >= K (deflation may result in N>K). \n
 * @param[in] N1
          N1 is INTEGER \n
          The location of the last eigenvalue in the leading submatrix.
          min(1,N) <= N1 <= N/2. \n
 * @param[out] D
          D is REAL array, dimension (N) \n
          D(I) contains the updated eigenvalues for
          1 <= I <= K. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ,N) \n
          Initially the first K columns are used as workspace.
          On output the columns 1 to K contain
          the updated eigenvectors. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= fla_max(1,N). \n
 * @param[in] RHO
          RHO is REAL \n
          The value of the parameter in the rank one update equation.
          RHO >= 0 required. \n
 * @param[in,out] DLAMDA
          DLAMDA is REAL array, dimension (K) \n
          The first K elements of this array contain the old roots
          of the deflated updating problem.  These are the poles
          of the secular equation. May be changed on output by
          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
          Cray-2, or Cray C-90, as described above. \n
 * @param[in] Q2
          Q2 is REAL array, dimension (LDQ2*N) \n
          The first K columns of this matrix contain the non-deflated
          eigenvectors for the split problem. \n
 * @param[in] INDX
          INDX is INTEGER array, dimension (N) \n
          The permutation used to arrange the columns of the deflated
          Q matrix into three groups (see SLAED2).
          The rows of the eigenvectors found by SLAED4 must be likewise
          permuted before the matrix multiply can take place. \n
 * @param[in] CTOT
          CTOT is INTEGER array, dimension (4) \n
          A count of the total number of the various types of columns
          in Q, as described in INDX.  The fourth column type is any
          column which has been deflated. \n
 * @param[in,out] W
          W is REAL array, dimension (K) \n
          The first K elements of this array contain the components
          of the deflation-adjusted updating vector. Destroyed on
          output. \n
 * @param[out] S
          S is REAL array, dimension (N1 + 1)*K \n
          Will contain the eigenvectors of the repaired matrix which
          will be multiplied by the previously accumulated eigenvectors
          to update the system. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, an eigenvalue did not converge \n

 * */
    template <typename T>
    void laed3(integer *k, integer *n, integer *n1, T *d, T *q, integer *ldq, T *rho, T *dlamda,
               T *q2, integer *indx, integer *ctot, T *w, T *s, integer *info)
    {
        laed3(*k, *n, *n1, d, q, *ldq, rho, dlamda, q2, *indx, *ctot, w, s, *info);
    }
    /** @}*/ // end of laed3

    /** @defgroup laed4 laed4
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED4 used by sstedc. Finds a single root of the secular equation

 * @details
 * \b Purpose:
    \verbatim
    This subroutine computes the I-th updated eigenvalue of a symmetric
    rank-one modification to a diagonal matrix whose elements are
    given in the array d, and that

               D(i) < D(j)  for  i < j

    and that RHO > 0.  This is arranged by the calling routine, and is
    no loss in generality.  The rank-one modified system is thus

               diag(D)  +  RHO * Z * Z_transpose.

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
          order, D(I) < D(J)  for I < J. \n
 * @param[in] Z
          Z is REAL array, dimension (N) \n
          The components of the updating vector. \n
 * @param[out] DELTA
          DELTA is REAL array, dimension (N) \n
          If N > 2, DELTA contains (D(j) - lambda_I) in its  j-th
          component.  If N = 1, then DELTA(1) = 1. If N = 2, see SLAED5
          for detail. The vector DELTA contains the information necessary
          to construct the eigenvectors by SLAED3 and SLAED9. \n
 * @param[in] RHO
          RHO is REAL \n
          The scalar in the symmetric updating formula. \n
 * @param[out] DLAM
          DLAM is REAL \n
          The computed lambda_I, the I-th updated eigenvalue. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          > 0:  if INFO = 1, the updating process failed.  \n

 * */
    template <typename T>
    void laed4(integer *n, integer *i, T *d, T *z, T *delta, T *rho, T *dlam, integer *info)
    {
        laed4(*n, *i, d, z, delta, rho, dlam, *info);
    }
    /** @}*/ // end of laed4

    /** @defgroup laed5 laed5
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED5 used by sstedc. Solves the 2-by-2 secular equation

 * @details
 * \b Purpose:
    \verbatim
    This subroutine computes the I-th eigenvalue of a symmetric rank-one
    modification of a 2-by-2 diagonal matrix

               diag(D)  +  RHO * Z * transpose(Z) .

    The diagonal elements in the array D are assumed to satisfy

               D(i) < D(j)  for  i < j .

    We also assume RHO > 0 and that the Euclidean norm of the vector
    Z is one.
    \endverbatim

 * @param[in] I
          I is INTEGER \n
          The index of the eigenvalue to be computed.  I = 1 or I = 2. \n
 * @param[in] D
          D is REAL array, dimension (2) \n
          The original eigenvalues.  We assume D(1) < D(2). \n
 * @param[in] Z
          Z is REAL array, dimension (2) \n
          The components of the updating vector. \n
 * @param[out] DELTA
          DELTA is REAL array, dimension (2) \n
          The vector DELTA contains the information necessary
          to construct the eigenvectors. \n
 * @param[in] RHO
          RHO is REAL \n
          The scalar in the symmetric updating formula. \n
 * @param[out] DLAM
          DLAM is REAL \n
          The computed lambda_I, the I-th updated eigenvalue. \n

 * */
    template <typename T>
    void laed5(integer *i, T *d, T *z, T *delta, T *rho, T *dlam)
    {
        laed5(*i, d, z, delta, rho, dlam);
    }
    /** @}*/ // end of laed5

    /** @defgroup laed6 laed6
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED6 used by sstedc. Computes one Newton step in solution of the secular
 equation

 * @details
 * \b Purpose:
    \verbatim
    LAED6 computes the positive or negative root (closest to the origin)
    of
                     z(1)        z(2)        z(3)
    f(x) =   rho + --------- + ---------- + ---------
                    d(1)-x      d(2)-x      d(3)-x

    It is assumed that

          if ORGATI = .true. the root is between d(2) and d(3);
          otherwise it is between d(1) and d(2)

    This routine will be called by SLAED4 when necessary. In most cases,
    the root sought is the smallest in magnitude, though it might not be
    in some extremely rare situations.
    \endverbatim

 * @param[in] KNITER
          KNITER is INTEGER \n
          Refer to SLAED4 for its significance. \n
 * @param[in] ORGATI
          ORGATI is LOGICAL \n
          If ORGATI is true, the needed root is between d(2) and
          d(3); otherwise it is between d(1) and d(2).  See
          SLAED4 for further details. \n
 * @param[in] RHO
          RHO is REAL \n
          Refer to the equation f(x) above. \n
 * @param[in] D
          D is REAL array, dimension (3) \n
          D satisfies d(1) < d(2) < d(3). \n
 * @param[in] Z
          Z is REAL array, dimension (3) \n
          Each of the elements in z must be positive. \n
 * @param[in] FINIT
          FINIT is REAL \n
          The value of f at 0. It is more accurate than the one
          evaluated inside this routine (if someone wants to do
          so). \n
 * @param[out] TAU
          TAU is REAL \n
          The root of the equation f(x). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          > 0: if INFO = 1, failure to converge \n

 * */
    template <typename T>
    void laed6(integer *kniter, logical *orgati, T *rho, T *d, T *z, T *finit, T *tau,
               integer *info)
    {
        laed6(kniter, orgati, rho, d, z, finit, tau, *info);
    }
    /** @}*/ // end of laed6

    /** @defgroup laed7 laed7
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED7 used by sstedc. Computes the updated eigensystem of a diagonal \n
     matrix after modification by a rank-one symmetric matrix. Used when the    \n
     original matrix is dense
 * @details
 * \b Purpose:
    \verbatim
    LAED7 computes the updated eigensystem of a diagonal
    matrix after modification by a rank-one symmetric matrix. This
    routine is used only for the eigenproblem which requires all
    eigenvalues and optionally eigenvectors of a dense symmetric matrix
    that has been reduced to tridiagonal form.  SLAED1 handles
    the case in which all eigenvalues and eigenvectors of a symmetric
    tridiagonal matrix are desired.

      T = Q(in) (D(in) + RHO * Z*Z**T) Q**T(in) = Q(out) * D(out) * Q**T(out)

       where Z = Q**Tu, u is a vector of length N with ones in the
       CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.

       The eigenvectors of the original matrix are stored in Q, and the
       eigenvalues are in D.  The algorithm consists of three stages:

          The first stage consists of deflating the size of the problem
          when there are multiple eigenvalues or if there is a zero in
          the Z vector.  For each such occurrence the dimension of the
          secular equation problem is reduced by one.  This stage is
          performed by the routine SLAED8.

          The second stage consists of calculating the updated
          eigenvalues. This is done by finding the roots of the secular
          equation via the routine SLAED4 (as called by SLAED9).
          This routine also calculates the eigenvectors of the current
          problem.

          The final stage consists of computing the updated eigenvectors
          directly using the updated eigenvalues.  The eigenvectors for
          the current problem are multiplied with the eigenvectors from
          the overall problem.
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          = 0:  Compute eigenvalues only. \n
          = 1:  Compute eigenvectors of original dense symmetric matrix
                also.  On entry, Q contains the orthogonal matrix used
                to reduce the original matrix to tridiagonal form. \n
 * @param[in] N
          N is INTEGER \n
          The dimension of the symmetric tridiagonal matrix.  N >= 0. \n
 * @param[in] QSIZ
          QSIZ is INTEGER \n
          The dimension of the orthogonal matrix used to reduce
          the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. \n
 * @param[in] TLVLS
          TLVLS is INTEGER \n
          The total number of merging levels in the overall divide and
          conquer tree. \n
 * @param[in] CURLVL
          CURLVL is INTEGER \n
          The current level in the overall merge routine,
          0 <= CURLVL <= TLVLS. \n
 * @param[in] CURPBM
          CURPBM is INTEGER \n
          The current problem in the current level in the overall
          merge routine (counting from upper left to lower right). \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the eigenvalues of the rank-1-perturbed matrix.
          On exit, the eigenvalues of the repaired matrix. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ, N) \n
          On entry, the eigenvectors of the rank-1-perturbed matrix.
          On exit, the eigenvectors of the repaired tridiagonal matrix. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= fla_max(1,N). \n
 * @param[out] INDXQ
          INDXQ is INTEGER array, dimension (N) \n
          The permutation which will reintegrate the subproblem just
          solved back into sorted order, i.e., D(INDXQ(I = 1, N))
          will be in ascending order. \n
 * @param[in] RHO
          RHO is REAL \n
          The subdiagonal element used to create the rank-1
          modification. \n
 * @param[in] CUTPNT
          CUTPNT is INTEGER \n
          Contains the location of the last eigenvalue in the leading
          sub-matrix.  min(1,N) <= CUTPNT <= N. \n
 * @param[in,out] QSTORE
          QSTORE is REAL array, dimension (N**2+1) \n
          Stores eigenvectors of submatrices encountered during
          divide and conquer, packed together. QPTR points to
          beginning of the submatrices. \n
 * @param[in,out] QPTR
          QPTR is INTEGER array, dimension (N+2) \n
          List of indices pointing to beginning of submatrices stored
          in QSTORE. The submatrices are numbered starting at the
          bottom left of the divide and conquer tree, from left to
          right and bottom to top. \n
 * @param[in] PRMPTR
          PRMPTR is INTEGER array, dimension (N lg N) \n
          Contains a list of pointers which indicate where in PERM a
          level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
          indicates the size of the permutation and also the size of
          the full, non-deflated problem. \n
 * @param[in] PERM
          PERM is INTEGER array, dimension (N lg N) \n
          Contains the permutations (from deflation and sorting) to be
          applied to each eigenblock. \n
 * @param[in] GIVPTR
          GIVPTR is INTEGER array, dimension (N lg N) \n
          Contains a list of pointers which indicate where in GIVCOL a
          level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
          indicates the number of Givens rotations. \n
 * @param[in] GIVCOL
          GIVCOL is INTEGER array, dimension (2, N lg N) \n
          Each pair of numbers indicates a pair of columns to take place
          in a Givens rotation. \n
 * @param[in] GIVNUM
          GIVNUM is REAL array, dimension (2, N lg N) \n
          Each number indicates the S value to be used in the
          corresponding Givens rotation. \n
 * @param[out] WORK
          WORK is REAL array, dimension (3*N+2*QSIZ*N) \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (4*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, an eigenvalue did not converge \n

 * */
    template <typename T>
    void laed7(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl,
               integer *curpbm, T *d, T *q, integer *ldq, integer *indxq, T *rho, integer *cutpnt,
               T *qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr,
               integer *givcol, T *givnum, T *work, integer *iwork, integer *info)
    {
        laed7(icompq, n, qsiz, tlvls, curlvl, curpbm, d, q, ldq, indxq, rho, cutpnt, qstore, qptr,
              prmptr, perm, givptr, givcol, givnum, work, iwork, info);
    }
    template <typename T, typename Ta>
    void laed7(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl,
               integer *curpbm, Ta *d, T *q, integer *ldq, Ta *rho, integer *indxq, Ta *qstore,
               integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol,
               Ta *givnum, T *work, Ta *rwork, integer *iwork, integer *info)
    {
        laed7(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr, prmptr,
              perm, givptr, givcol, givnum, work, rwork, iwork, info);
    }
    /** @}*/ // end of laed7

    /** @defgroup laed8 laed8
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED8 used by sstedc. Merges eigenvalues and deflates secular \n
     equation. Used when the original matrix is dense
 * @details
 * \b Purpose:
    \verbatim
    LAED8 merges the two sets of eigenvalues together into a single
    sorted set.  Then it tries to deflate the size of the problem.
    There are two ways in which deflation can occur:  when two or more
    eigenvalues are close together or if there is a tiny element in the
    Z vector.  For each such occurrence the order of the related secular
    equation problem is reduced by one.
     \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          = 0:  Compute eigenvalues only. \n
          = 1:  Compute eigenvectors of original dense symmetric matrix
                also.  On entry, Q contains the orthogonal matrix used
                to reduce the original matrix to tridiagonal form. \n
 * @param[out] K
          K is INTEGER \n
          The number of non-deflated eigenvalues, and the order of the
          related secular equation. \n
 * @param[in] N
          N is INTEGER \n
          The dimension of the symmetric tridiagonal matrix.  N >= 0. \n
 * @param[in] QSIZ
          QSIZ is INTEGER \n
          The dimension of the orthogonal matrix used to reduce
          the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the eigenvalues of the two submatrices to be
          combined.  On exit, the trailing (N-K) updated eigenvalues
          (those which were deflated) sorted into increasing order. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          If ICOMPQ = 0, Q is not referenced.  Otherwise,
          on entry, Q contains the eigenvectors of the partially solved
          system which has been previously updated in matrix
          multiplies with other partially solved eigensystems.
          On exit, Q contains the trailing (N-K) updated eigenvectors
          (those which were deflated) in its last N-K columns. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= fla_max(1,N). \n
 * @param[in] INDXQ
          INDXQ is INTEGER array, dimension (N) \n
          The permutation which separately sorts the two sub-problems
          in D into ascending order.  Note that elements in the second
          half of this permutation must first have CUTPNT added to
          their values in order to be accurate. \n
 * @param[in,out] RHO
          RHO is REAL \n
          On entry, the off-diagonal element associated with the rank-1
          cut which originally split the two submatrices which are now
          being recombined. \n
          On exit, RHO has been modified to the value required by
          SLAED3. \n
 * @param[in] CUTPNT
          CUTPNT is INTEGER \n
          The location of the last eigenvalue in the leading
          sub-matrix.  min(1,N) <= CUTPNT <= N. \n
 * @param[in] Z
          Z is REAL array, dimension (N) \n
          On entry, Z contains the updating vector (the last row of
          the first sub-eigenvector matrix and the first row of the
          second sub-eigenvector matrix). \n
          On exit, the contents of Z are destroyed by the updating
          process. \n
 * @param[out] DLAMDA
          DLAMDA is REAL array, dimension (N) \n
          A copy of the first K eigenvalues which will be used by
          SLAED3 to form the secular equation. \n
 * @param[out] Q2
          Q2 is REAL array, dimension (LDQ2,N) \n
          If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
          a copy of the first K eigenvectors which will be used by
          SLAED7 in a matrix multiply (SGEMM) to update the new
          eigenvectors. \n
 * @param[in] LDQ2
          LDQ2 is INTEGER \n
          The leading dimension of the array Q2.  LDQ2 >= fla_max(1,N). \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first k values of the final deflation-altered z-vector and
          will be passed to SLAED3. \n
 * @param[out] PERM
          PERM is INTEGER array, dimension (N) \n
          The permutations (from deflation and sorting) to be applied
          to each eigenblock. \n
 * @param[out] GIVPTR
          GIVPTR is INTEGER \n
          The number of Givens rotations which took place in this
          subproblem. \n
 * @param[out] GIVCOL
          GIVCOL is INTEGER array, dimension (2, N) \n
          Each pair of numbers indicates a pair of columns to take place
          in a Givens rotation. \n
 * @param[out] GIVNUM
          GIVNUM is REAL array, dimension (2, N) \n
          Each number indicates the S value to be used in the
          corresponding Givens rotation. \n
 * @param[out] INDXP
          INDXP is INTEGER array, dimension (N) \n
          The permutation used to place deflated values of D at the end
          of the array.  INDXP(1:K) points to the nondeflated D-values
          and INDXP(K+1:N) points to the deflated eigenvalues. \n
 * @param[out] INDX
          INDX is INTEGER array, dimension (N) \n
          The permutation used to sort the contents of D into ascending
          order. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 * */
    template <typename T>
    void laed8(integer *icompq, integer *k, integer *n, integer *qsiz, T *d, T *q, integer *ldq,
               integer *indxq, T *rho, integer *cutpnt, T *z, T *dlamda, T *q2, integer *ldq2, T *w,
               integer *perm, integer *givptr, integer *givcol, T *givnum, integer *indxp,
               integer *indx, integer *info)
    {
        laed8(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, z, dlamda, q2, ldq2, w, perm,
              givptr, givcol, givnum, indxp, indx, info);
    }
    template <typename T, typename Ta>
    void laed8(integer *k, integer *n, integer *qsiz, T *q, integer *ldq, Ta *d, Ta *rho,
               integer *cutpnt, Ta *z, Ta *dlamda, T *q2, integer *ldq2, Ta *w, integer *indxp,
               integer *indx, integer *indxq, integer *perm, integer *givptr, integer *givcol,
               Ta *givnum, integer *info)
    {
        laed8(k, n, qsiz, q, ldq, d, rho, cutpnt, z, dlamda, q2, ldq2, w, indxp, indx, indxq, perm,
              givptr, givcol, givnum, info);
    }
    /** @}*/ // end of laed8

    /** @defgroup laed9 laed9
     * @ingroup SYM
     * @{
     */
    /*! @brief LAED9 used by sstedc. Finds the roots of the secular equation and \n
     updates the eigenvectors. Used when the original matrix is dense
 * @details
 * \b Purpose:
    \verbatim
    LAED9 finds the roots of the secular equation, as defined by the
    values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
    appropriate calls to SLAED4 and then stores the new matrix of
    eigenvectors for use in calculating the next level of Z vectors.
    \endverbatim

 * @param[in] K
          K is INTEGER \n
          The number of terms in the rational function to be solved by
          SLAED4.  K >= 0. \n
 * @param[in] KSTART
          KSTART is INTEGER \n
 * @param[in] KSTOP
          KSTOP is INTEGER \n
          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP
          are to be computed.  1 <= KSTART <= KSTOP <= K. \n
 * @param[in] N
          N is INTEGER \n
          The number of rows and columns in the Q matrix.
          N >= K (delation may result in N > K). \n
 * @param[out] D
          D is REAL array, dimension (N) \n
          D(I) contains the updated eigenvalues \n
          for KSTART <= I <= KSTOP. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ,N) \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= fla_max(1, N). \n
 * @param[in] RHO
          RHO is REAL \n
          The value of the parameter in the rank one update equation.
          RHO >= 0 required. \n
 * @param[in] DLAMDA
          DLAMDA is REAL array, dimension (K) \n
          The first K elements of this array contain the old roots
          of the deflated updating problem.  These are the poles
          of the secular equation. \n
 * @param[in] W
          W is REAL array, dimension (K) \n
          The first K elements of this array contain the components
          of the deflation-adjusted updating vector. \n
 * @param[out] S
          S is REAL array, dimension (LDS, K) \n
          Will contain the eigenvectors of the repaired matrix which
          will be stored for subsequent Z vector calculation and
          multiplied by the previously accumulated eigenvectors
          to update the system. \n
 * @param[in] LDS
          LDS is INTEGER \n
          The leading dimension of S.  LDS >= fla_max(1, K). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, an eigenvalue did not converge \n

 * */
    template <typename T>
    void laed9(integer *k, integer *kstart, integer *kstop, integer *n, T *d, T *q, integer *ldq,
               T *rho, T *dlamda, T *w, T *s, integer *lds, integer *info)
    {
        laed9(k, kstart, kstop, n, d, q, ldq, rho, dlamda, w, s, lds, info);
    }
    /** @}*/ // end of laed9

    /** @defgroup laeda laeda
     * @ingroup SYM
     * @{
     */
    /*! @brief LAEDA used by sstedc. Computes the Z vector determining the rank-one \n
     modification of the diagonal matrix. Used when the original matrix is dense
 * @details
 * \b Purpose:
    \verbatim
    LAEDA computes the Z vector corresponding to the merge step in the
    CURLVLth step of the merge process with TLVLS steps for the CURPBMth
    problem.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The dimension of the symmetric tridiagonal matrix.  N >= 0. \n
 * @param[in] TLVLS
          TLVLS is INTEGER \n
          The total number of merging levels in the overall divide and
          conquer tree. \n
 * @param[in] CURLVL
          CURLVL is INTEGER \n
          The current level in the overall merge routine,
          0 <= curlvl <= tlvls. \n
 * @param[in] CURPBM
          CURPBM is INTEGER \n
          The current problem in the current level in the overall
          merge routine (counting from upper left to lower right). \n
 * @param[in] PRMPTR
          PRMPTR is INTEGER array, dimension (N lg N) \n
          Contains a list of pointers which indicate where in PERM a
          level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
          indicates the size of the permutation and incidentally the
          size of the full, non-deflated problem. \n
 * @param[in] PERM
          PERM is INTEGER array, dimension (N lg N) \n
          Contains the permutations (from deflation and sorting) to be
          applied to each eigenblock. \n
 * @param[in] GIVPTR
          GIVPTR is INTEGER array, dimension (N lg N) \n
          Contains a list of pointers which indicate where in GIVCOL a
          level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
          indicates the number of Givens rotations. \n
 * @param[in] GIVCOL
          GIVCOL is INTEGER array, dimension (2, N lg N) \n
          Each pair of numbers indicates a pair of columns to take place
          in a Givens rotation. \n
 * @param[in] GIVNUM
          GIVNUM is REAL array, dimension (2, N lg N) \n
          Each number indicates the S value to be used in the
          corresponding Givens rotation. \n
 * @param[in] Q
          Q is REAL array, dimension (N**2) \n
          Contains the square eigenblocks from previous levels, the
          starting positions for blocks are given by QPTR. \n
 * @param[in] QPTR
          QPTR is INTEGER array, dimension (N+2) \n
          Contains a list of pointers which indicate where in Q an
          eigenblock is stored.  SQRT(QPTR(i+1) - QPTR(i)) indicates
          the size of the block. \n
 * @param[out] Z
          Z is REAL array, dimension (N) \n
          On output this vector contains the updating vector (the last
          row of the first sub-eigenvector matrix and the first row of
          the second sub-eigenvector matrix). \n
 * @param[out] ZTEMP
          ZTEMP is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 * */
    template <typename T>
    void laeda(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr,
               integer *perm, integer *givptr, integer *givcol, float *givnum, float *q,
               integer *qptr, float *z, float *ztemp, integer *info)
    {
        laeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp,
              info);
    }
    /** @}*/ // end of laeda

    /** @defgroup larra larra
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRA computes the splitting points with the specified threshold

 * @details
 * \b Purpose:
    \verbatim
    Compute the splitting points with threshold SPLTOL.
    LARRA sets any "small" off-diagonal elements to zero.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix. N > 0. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          On entry, the N diagonal elements of the tridiagonal
          matrix T. \n
 * @param[in,out] E
          E is REAL array, dimension (N) \n
          On entry, the first (N-1) entries contain the subdiagonal
          elements of the tridiagonal matrix T; E(N) need not be set.
          On exit, the entries E(ISPLIT(I)), 1 <= I <= NSPLIT,
          are set to zero, the other entries of E are untouched. \n
 * @param[in,out] E2
          E2 is REAL array, dimension (N) \n
          On entry, the first (N-1) entries contain the SQUARES of the
          subdiagonal elements of the tridiagonal matrix T;
          E2(N) need not be set. \n
          On exit, the entries E2(ISPLIT(I)),
          1 <= I <= NSPLIT, have been set to zero \n
 * @param[in] SPLTOL
          SPLTOL is REAL \n
          The threshold for splitting. Two criteria can be used:
          SPLTOL<0 : criterion based on absolute off-diagonal value
          SPLTOL>0 : criterion that preserves relative accuracy \n
 * @param[in] TNRM
          TNRM is REAL \n
          The norm of the matrix. \n
 * @param[out] NSPLIT
          NSPLIT is INTEGER \n
          The number of blocks T splits into. 1 <= NSPLIT <= N. \n
 * @param[out] ISPLIT
          ISPLIT is INTEGER array, dimension (N) \n
          The splitting points, at which T breaks up into blocks.
          The first block consists of rows/columns 1 to ISPLIT(1),
          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
          etc., and the NSPLIT-th consists of rows/columns
          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n

 *  * */
    template <typename T>
    void larra(integer *n, T *d, T *e, T *e2, T *spltol, T *tnrm, integer *nsplit, integer *isplit,
               integer *info)
    {
        larra(n, d, e, e2, spltol, tnrm, nsplit, isplit, info);
    }
    /** @}*/ // end of larra

    /** @defgroup larrb larrb
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRB provides limited bisection to locate eigenvalues for more accuracy

 * @details
 * \b Purpose:
    \verbatim
    Given the relatively robust representation(RRR) L D L^T, SLARRB
    does "limited" bisection to refine the eigenvalues of L D L^T,
    W(IFIRST-OFFSET) through W(ILAST-OFFSET), to more accuracy. Initial
    guesses for these eigenvalues are input in W, the corresponding estimate
    of the error in these guesses and their gaps are input in WERR
    and WGAP, respectively. During bisection, intervals
    [left, right] are maintained by storing their mid-points and
    semi-widths in the arrays W and WERR respectively.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The N diagonal elements of the diagonal matrix D. \n
 * @param[in] LLD
          LLD is REAL array, dimension (N-1) \n
          The (N-1) elements L(i)*L(i)*D(i). \n
 * @param[in] IFIRST
          IFIRST is INTEGER \n
          The index of the first eigenvalue to be computed. \n
 * @param[in] ILAST
          ILAST is INTEGER \n
          The index of the last eigenvalue to be computed. \n
 * @param[in] RTOL1
          RTOL1 is REAL \n
 * @param[in] RTOL2
          RTOL2 is REAL \n
          Tolerance for the convergence of the bisection intervals.
          An interval [LEFT,RIGHT] has converged if
          RIGHT-LEFT < MAX(RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|))
          where GAP is the (estimated) distance to the nearest
          eigenvalue. \n
 * @param[in] OFFSET
          OFFSET is INTEGER \n
          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET
          through ILAST-OFFSET elements of these arrays are to be used. \n
 * @param[in,out] W
          W is REAL array, dimension (N) \n
          On input, W(IFIRST-OFFSET) through W(ILAST-OFFSET) are
          estimates of the eigenvalues of L D L^T indexed IFIRST through
          ILAST. \n
          On output, these estimates are refined. \n
 * @param[in,out] WGAP
          WGAP is REAL array, dimension (N-1) \n
          On input, the (estimated) gaps between consecutive
          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between
          eigenvalues I and I+1. Note that if IFIRST = ILAST
          then WGAP(IFIRST-OFFSET) must be set to ZERO. \n
          On output, these gaps are refined. \n
 * @param[in,out] WERR
          WERR is REAL array, dimension (N) \n
          On input, WERR(IFIRST-OFFSET) through WERR(ILAST-OFFSET) are
          the errors in the estimates of the corresponding elements in W.
          On output, these errors are refined. \n
 * @param[out] WORK
          WORK is REAL array, dimension (2*N) \n
          Workspace. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (2*N) \n
          Workspace. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot in the Sturm sequence. \n
 * @param[in] SPDIAM
          SPDIAM is REAL \n
          The spectral diameter of the matrix. \n
 * @param[in] TWIST
          TWIST is INTEGER \n
          The twist index for the twisted factorization that is used
          for the negcount. \n
          TWIST = N: Compute negcount from L D L^T - LAMBDA I = L+ D+ L+^T \n
          TWIST = 1: Compute negcount from L D L^T - LAMBDA I = U- D- U-^T \n
          TWIST = R: Compute negcount from L D L^T - LAMBDA I = N(r) D(r) N(r) \n
 * @param[out] INFO
          INFO is INTEGER \n
          Error flag.  \n

 *  * */
    template <typename T>
    void larrb(integer *n, T *d, T *lld, integer *ifirst, integer *ilast, T *rtol1, T *rtol2,
               integer *offset, T *w, T *wgap, T *werr, T *work, integer *iwork, T *pivmin,
               T *spdiam, integer *twist, integer *info)
    {
        larrb(n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork, pivmin,
              spdiam, twist, info);
    }
    /** @}*/ // end of larrb

    /** @defgroup larrc larrc
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRC computes the number of eigenvalues of the symmetric tridiagonal matrix

 * @details
 * \b Purpose:
    \verbatim
    Find the number of eigenvalues of the symmetric tridiagonal matrix T
    that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T
    if JOBT = 'L'.
    \endverbatim

 * @param[in] JOBT
          JOBT is CHARACTER*1 \n
          = 'T':  Compute Sturm count for matrix T. \n
          = 'L':  Compute Sturm count for matrix L D L^T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix. N > 0. \n
 * @param[in] VL
          VL is REAL \n
          The lower bound for the eigenvalues. \n
 * @param[in] VU
          VU is REAL \n
          The upper bound for the eigenvalues. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T. \n
          JOBT = 'L': The N diagonal elements of the diagonal matrix D. \n
 * @param[in] E
          E is REAL array, dimension (N) \n
          JOBT = 'T': The N-1 offdiagonal elements of the matrix T. \n
          JOBT = 'L': The N-1 offdiagonal elements of the matrix L. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot in the Sturm sequence for T. \n
 * @param[out] EIGCNT
          EIGCNT is INTEGER \n
          The number of eigenvalues of the symmetric tridiagonal matrix T
          that are in the interval (VL,VU] \n
 * @param[out] LCNT
          LCNT is INTEGER \n
 * @param[out] RCNT
          RCNT is INTEGER \n
          The left and right negcounts of the interval. \n
 * @param[out] INFO
          INFO is INTEGER  \n

 *  * */
    template <typename T>
    void larrc(char *jobt, integer *n, T *vl, T *vu, T *d, T *e, T *pivmin, integer *eigcnt,
               integer *lcnt, integer *rcnt, integer *info)
    {
        larrc(jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info);
    }
    /** @}*/ // end of larrc

    /** @defgroup larrd larrd
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRD computes the eigenvalues of a symmetric tridiagonal matrix to suitable
 accuracy

 * @details
 * \b Purpose:
    \verbatim
    LARRD computes the eigenvalues of a symmetric tridiagonal
    matrix T to suitable accuracy. This is an auxiliary code to be
    called from SSTEMR.
    The user may ask for all eigenvalues, all eigenvalues
    in the half-open interval (VL, VU], or the IL-th through IU-th
    eigenvalues.

    To avoid overflow, the matrix must be scaled so that its
    largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and
 for greatest accuracy, it should not be much smaller than that.

    See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
    Matrix", Report CS41, Computer Science Dept., Stanford
    University, July 21, 1966.
    \endverbatim

 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': ("All")   all eigenvalues will be found. \n
          = 'V': ("Value") all eigenvalues in the half-open interval
                           (VL, VU] will be found. \n
          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
                           entire matrix) will be found. \n
 * @param[in] ORDER
          ORDER is CHARACTER*1 \n
          = 'B': ("By Block") the eigenvalues will be grouped by
                              split-off block (see IBLOCK, ISPLIT) and
                              ordered from smallest to largest within
                              the block. \n
          = 'E': ("Entire matrix") \n
                              the eigenvalues for the entire matrix
                              will be ordered from smallest to
                              largest. \n
 * @param[in] N
          N is INTEGER \n
          The order of the tridiagonal matrix T.  N >= 0. \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues.  Eigenvalues less than or equal
          to VL, or greater than VU, will not be   returned.  VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues.  Eigenvalues less than or equal
          to VL, or greater than VU, will not be   returned.  VL < VU. \n
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be   returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be   returned. \n
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] GERS
          GERS is REAL array, dimension (2*N) \n
          The N Gerschgorin intervals (the i-th Gerschgorin interval
          is (GERS(2*i-1), GERS(2*i)). \n
 * @param[in] RELTOL
          RELTOL is REAL \n
          The minimum relative width of an interval.  When an interval
          is narrower than RELTOL times the larger (in
          magnitude) endpoint, then it is considered to be
          sufficiently small, i.e., converged.  Note: this should
          always be at least radix*machine epsilon. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the tridiagonal matrix T. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) off-diagonal elements of the tridiagonal matrix T. \n
 * @param[in] E2
          E2 is REAL array, dimension (N-1) \n
          The (n-1) squared off-diagonal elements of the tridiagonal matrix T. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot allowed in the Sturm sequence for T. \n
 * @param[in] NSPLIT
          NSPLIT is INTEGER \n
          The number of diagonal blocks in the matrix T.
          1 <= NSPLIT <= N. \n
 * @param[in] ISPLIT
          ISPLIT is INTEGER array, dimension (N) \n
          The splitting points, at which T breaks up into submatrices.
          The first submatrix consists of rows/columns 1 to ISPLIT(1),
          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
          etc., and the NSPLIT-th consists of rows/columns
          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. \n
          (Only the first NSPLIT elements will actually be used, but
          since the user cannot know a priori what value NSPLIT will
          have, N words must be reserved for ISPLIT.) \n
 * @param[out] M
          M is INTEGER \n
          The actual number of eigenvalues found. 0 <= M <= N.
          (See also the description of INFO=2,3.) \n
 * @param[out] W
          W is REAL array, dimension (N) \n
          On exit, the first M elements of W will contain the
          eigenvalue approximations. SLARRD computes an interval
          I_j = (a_j, b_j] that includes eigenvalue j. The eigenvalue
          approximation is given as the interval midpoint
          W(j)= (a_j + b_j)/2. The corresponding error is bounded by
          WERR(j) = abs(a_j - b_j)/2 \n
 * @param[out] WERR
          WERR is REAL array, dimension (N) \n
          The error bound on the corresponding eigenvalue approximation
          in W. \n
 * @param[out] WL
          WL is REAL \n
 * @param[out] WU
          WU is REAL \n
          The interval (WL, WU] contains all the wanted eigenvalues. \n
          If RANGE='V', then WL=VL and WU=VU. \n
          If RANGE='A', then WL and WU are the global Gerschgorin bounds
                        on the spectrum. \n
          If RANGE='I', then WL and WU are computed by SLAEBZ from the
                        index range specified. \n
 * @param[out] IBLOCK
          IBLOCK is INTEGER array, dimension (N) \n
          At each row/column j where E(j) is zero or small, the
          matrix T is considered to split into a block diagonal
          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
          block (from 1 to the number of blocks) the eigenvalue W(i)
          belongs.  (SLARRD may use the remaining N-M elements as
          workspace.) \n
 * @param[out] INDEXW
          INDEXW is INTEGER array, dimension (N) \n
          The indices of the eigenvalues within each block (submatrix);
          for example, INDEXW(i)= j and IBLOCK(i)=k imply that the
          i-th eigenvalue W(i) is the j-th eigenvalue in block k. \n
 * @param[out] WORK
          WORK is REAL array, dimension (4*N) \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (3*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  some or all of the eigenvalues failed to converge or
                were not computed: \n
                =1 or 3: Bisection failed to converge for some
                        eigenvalues; these eigenvalues are flagged by a
                        negative block number.  The effect is that the
                        eigenvalues may not be as accurate as the
                        absolute and relative tolerances.  This is
                        generally caused by unexpectedly inaccurate
                        arithmetic. \n
                =2 or 3: RANGE='I' only: Not all of the eigenvalues
                        IL:IU were found. \n
                        Effect: M < IU+1-IL \n
                        Cause:  non-monotonic arithmetic, causing the
                                Sturm sequence to be non-monotonic. \n
                        Cure:   recalculate, using RANGE='A', and pick
                                out eigenvalues IL:IU.  In some cases,
                                increasing the PARAMETER "FUDGE" may
                                make things work. \n
                = 4:    RANGE='I', and the Gershgorin interval
                        initially used was too small.  No eigenvalues
                        were computed. \n
                        Probable cause: your machine has sloppy
                                        floating-point arithmetic. \n
                        Cure: Increase the PARAMETER "FUDGE",
                              recompile, and try again.  \n

 *  * */
    template <typename T>
    void larrd(char *range, char *order, integer *n, T *vl, T *vu, integer *il, integer *iu,
               T *gers, T *reltol, T *d, T *e, T *e2, T *pivmin, integer *nsplit, integer *isplit,
               integer *m, T *w, T *werr, T *wl, T *wu, integer *iblock, integer *indexw, T *work,
               integer *iwork, integer *info)
    {
        larrd(range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w,
              werr, wl, wu, iblock, indexw, work, iwork, info);
    }
    /** @}*/ // end of larrd

    /** @defgroup larre larre
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRE given the tridiagonal matrix T, sets small off-diagonal elements to \n
     zero and for each unreduced block Ti, finds base representations and eigenvalues
 * @details
 * \b Purpose:
    \verbatim
    To find the desired eigenvalues of a given real symmetric
    tridiagonal matrix T, SLARRE sets any "small" off-diagonal
    elements to zero, and for each unreduced block T_i, it finds
    (a) a suitable shift at one end of the block's spectrum,
    (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and
    (c) eigenvalues of each L_i D_i L_i^T.
    The representations and eigenvalues found are then used by
    SSTEMR to compute the eigenvectors of T.
    The accuracy varies depending on whether bisection is used to
    find a few eigenvalues or the dqds algorithm (subroutine SLASQ2) to
    conpute all and then discard any unwanted one.
    As an added benefit, SLARRE also outputs the n
    Gerschgorin intervals for the matrices L_i D_i L_i^T.
    \endverbatim

 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': ("All")   all eigenvalues will be found. \n
          = 'V': ("Value") all eigenvalues in the half-open interval
                           (VL, VU] will be found. \n
          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
                           entire matrix) will be found. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix. N > 0. \n
 * @param[in,out] VL
          VL is REAL \n
          If RANGE='V', the lower bound for the eigenvalues.
          Eigenvalues less than or equal to VL, or greater than VU,
          will not be   returned.  VL < VU. \n
          If RANGE='I' or ='A', SLARRE computes bounds on the desired
          part of the spectrum. \n
 * @param[in,out] VU
          VU is REAL \n
          If RANGE='V', the upper bound for the eigenvalues.
          Eigenvalues less than or equal to VL, or greater than VU,
          will not be   returned.  VL < VU. \n
          If RANGE='I' or ='A', SLARRE computes bounds on the desired
          part of the spectrum. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest eigenvalue to be   returned. \n
          1 <= IL <= IU <= N. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest eigenvalue to be   returned. \n
          1 <= IL <= IU <= N. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the N diagonal elements of the tridiagonal
          matrix T. \n
          On exit, the N diagonal elements of the diagonal
          matrices D_i. \n
 * @param[in,out] E
          E is REAL array, dimension (N) \n
          On entry, the first (N-1) entries contain the subdiagonal
          elements of the tridiagonal matrix T; E(N) need not be set. \n
          On exit, E contains the subdiagonal elements of the unit
          bidiagonal matrices L_i. The entries E(ISPLIT(I)), \n
          1 <= I <= NSPLIT, contain the base points sigma_i on output. \n
 * @param[in,out] E2
          E2 is REAL array, dimension (N) \n
          On entry, the first (N-1) entries contain the SQUARES of the
          subdiagonal elements of the tridiagonal matrix T;
          E2(N) need not be set. \n
          On exit, the entries E2(ISPLIT(I)),
          1 <= I <= NSPLIT, have been set to zero \n
 * @param[in] RTOL1
          RTOL1 is REAL \n
 * @param[in] RTOL2
          RTOL2 is REAL \n
          Parameters for bisection. \n
          An interval [LEFT,RIGHT] has converged if
          RIGHT-LEFT < MAX(RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|)) \n
 * @param[in] SPLTOL
          SPLTOL is REAL \n
          The threshold for splitting. \n
 * @param[out] NSPLIT
          NSPLIT is INTEGER \n
          The number of blocks T splits into. 1 <= NSPLIT <= N. \n
 * @param[out] ISPLIT
          ISPLIT is INTEGER array, dimension (N) \n
          The splitting points, at which T breaks up into blocks.
          The first block consists of rows/columns 1 to ISPLIT(1),
          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
          etc., and the NSPLIT-th consists of rows/columns
          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. \n
 * @param[out] M
          M is INTEGER \n
          The total number of eigenvalues (of all L_i D_i L_i^T) \n
          found.
 * @param[out] W
          W is REAL array, dimension (N) \n
          The first M elements contain the eigenvalues. The
          eigenvalues of each of the blocks, L_i D_i L_i^T, are
          sorted in ascending order (SLARRE may use the
          remaining N-M elements as workspace). \n
 * @param[out] WERR
          WERR is REAL array, dimension (N) \n
          The error bound on the corresponding eigenvalue in W. \n
 * @param[out] WGAP
          WGAP is REAL array, dimension (N) \n
          The separation from the right neighbor eigenvalue in W.
          The gap is only with respect to the eigenvalues of the same block
          as each block has its own representation tree.
          Exception: at the right end of a block we store the left gap \n
 * @param[out] IBLOCK
          IBLOCK is INTEGER array, dimension (N) \n
          The indices of the blocks (submatrices) associated with the
          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue
          W(i) belongs to the first block from the top, =2 if W(i)
          belongs to the second block, etc. \n
 * @param[out] INDEXW
          INDEXW is INTEGER array, dimension (N) \n
          The indices of the eigenvalues within each block (submatrix);
          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the
          i-th eigenvalue W(i) is the 10-th eigenvalue in block 2 \n
 * @param[out] GERS
          GERS is REAL array, dimension (2*N) \n
          The N Gerschgorin intervals (the i-th Gerschgorin interval
          is (GERS(2*i-1), GERS(2*i)). \n
 * @param[out] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot in the Sturm sequence for T. \n
 * @param[out] WORK
          WORK is REAL array, dimension (6*N) \n
          Workspace. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (5*N) \n
          Workspace. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          > 0:  A problem occurred in SLARRE. \n
          < 0:  One of the called subroutines signaled an internal problem.
                Needs inspection of the corresponding parameter IINFO
                for further information. \n
          =-1:  Problem in SLARRD. \n
          = 2:  No base representation could be found in MAXTRY iterations.
                Increasing MAXTRY and recompilation might be a remedy. \n
          =-3:  Problem in SLARRB when computing the refined root
                representation for SLASQ2. \n
          =-4:  Problem in SLARRB when preforming bisection on the
                desired part of the spectrum. \n
          =-5:  Problem in SLASQ2. \n
          =-6:  Problem in SLASQ2.  \n

 *  * */
    template <typename T>
    void larre(char *range, integer *n, T *vl, T *vu, integer *il, integer *iu, T *d, T *e, T *e2,
               T *rtol1, T *rtol2, T *spltol, integer *nsplit, integer *isplit, integer *m, T *w,
               T *werr, T *wgap, integer *iblock, integer *indexw, T *gers, T *pivmin, T *work,
               integer *iwork, integer *info)
    {
        larre(range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr,
              wgap, iblock, indexw, gers, pivmin, work, iwork, info);
    }
    /** @}*/ // end of larre

    /** @defgroup larrf larrf
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRF finds a new relatively robust representation such that at least \n
     one of the eigenvalues is relatively isolated
 * @details
 * \b Purpose:
    \verbatim
    Given the initial representation L D L^T and its cluster of close
    eigenvalues (in a relative measure), W(CLSTRT), W(CLSTRT+1), ...
    W(CLEND), SLARRF finds a new relatively robust representation
    L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the
    eigenvalues of L(+) D(+) L(+)^T is relatively isolated.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix (subblock, if the matrix split). \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The N diagonal elements of the diagonal matrix D. \n
 * @param[in] L
          L is REAL array, dimension (N-1) \n
          The (N-1) subdiagonal elements of the unit bidiagonal
          matrix L. \n
 * @param[in] LD
          LD is REAL array, dimension (N-1) \n
          The (N-1) elements L(i)*D(i). \n
 * @param[in] CLSTRT
          CLSTRT is INTEGER \n
          The index of the first eigenvalue in the cluster. \n
 * @param[in] CLEND
          CLEND is INTEGER \n
          The index of the last eigenvalue in the cluster. \n
 * @param[in] W
          W is REAL array, dimension \n
          dimension is >=  (CLEND-CLSTRT+1) \n
          The eigenvalue APPROXIMATIONS of L D L^T in ascending order.
          W(CLSTRT) through W(CLEND) form the cluster of relatively
          close eigenalues. \n
 * @param[in,out] WGAP
          WGAP is REAL array, dimension \n
          dimension is >=  (CLEND-CLSTRT+1) \n
          The separation from the right neighbor eigenvalue in W. \n
 * @param[in] WERR
          WERR is REAL array, dimension \n
          dimension is >=  (CLEND-CLSTRT+1) \n
          WERR contain the semiwidth of the uncertainty
          interval of the corresponding eigenvalue APPROXIMATION in W \n
 * @param[in] SPDIAM
          SPDIAM is REAL \n
          estimate of the spectral diameter obtained from the
          Gerschgorin intervals \n
 * @param[in] CLGAPL
          CLGAPL is REAL \n
 * @param[in] CLGAPR
          CLGAPR is REAL \n
          absolute gap on each end of the cluster.
          Set by the calling routine to protect against shifts too close
          to eigenvalues outside the cluster. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot allowed in the Sturm sequence. \n
 * @param[out] SIGMA
          SIGMA is REAL \n
          The shift used to form L(+) D(+) L(+)^T. \n
 * @param[out] DPLUS
          DPLUS is REAL array, dimension (N) \n
          The N diagonal elements of the diagonal matrix D(+). \n
 * @param[out] LPLUS
          LPLUS is REAL array, dimension (N-1) \n
          The first (N-1) elements of LPLUS contain the subdiagonal
          elements of the unit bidiagonal matrix L(+). \n
 * @param[out] WORK
          WORK is REAL array, dimension (2*N) \n
          Workspace. \n
 * @param[out] INFO
          INFO is INTEGER \n
          Signals processing OK (=0) or failure (=1) \n

 *  * */
    template <typename T>
    void larrf(integer *n, T *d, T *l, T *ld, integer *clstrt, integer *clend, T *w, T *wgap,
               T *werr, T *spdiam, T *clgapl, T *clgapr, T *pivmin, T *sigma, T *dplus, T *lplus,
               T *work, integer *info)
    {
        larrf(n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, sigma,
              dplus, lplus, work, info);
    }
    /** @}*/ // end of larrf

    /** @defgroup larrj larrj
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRJ performs refinement of the initial estimates of the eigenvalues of the
 matrix T

 * @details
 * \b Purpose:
    \verbatim
    Given the initial eigenvalue approximations of T, SLARRJ
    does  bisection to refine the eigenvalues of T,
    W(IFIRST-OFFSET) through W(ILAST-OFFSET), to more accuracy. Initial
    guesses for these eigenvalues are input in W, the corresponding estimate
    of the error in these guesses in WERR. During bisection, intervals
    [left, right] are maintained by storing their mid-points and
    semi-widths in the arrays W and WERR respectively.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The N diagonal elements of T. \n
 * @param[in] E2
          E2 is REAL array, dimension (N-1) \n
          The Squares of the (N-1) subdiagonal elements of T. \n
 * @param[in] IFIRST
          IFIRST is INTEGER \n
          The index of the first eigenvalue to be computed. \n
 * @param[in] ILAST
          ILAST is INTEGER \n
          The index of the last eigenvalue to be computed. \n
 * @param[in] RTOL
          RTOL is REAL \n
          Tolerance for the convergence of the bisection intervals.
          An interval [LEFT,RIGHT] has converged if
          RIGHT-LEFT < RTOL*MAX(|LEFT|,|RIGHT|). \n
 * @param[in] OFFSET
          OFFSET is INTEGER \n
          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET
          through ILAST-OFFSET elements of these arrays are to be used. \n
 * @param[in,out] W
          W is REAL array, dimension (N) \n
          On input, W(IFIRST-OFFSET) through W(ILAST-OFFSET) are
          estimates of the eigenvalues of L D L^T indexed IFIRST through
          ILAST. \n
          On output, these estimates are refined. \n
 * @param[in,out] WERR
          WERR is REAL array, dimension (N) \n
          On input, WERR(IFIRST-OFFSET) through WERR(ILAST-OFFSET) are
          the errors in the estimates of the corresponding elements in W.
          On output, these errors are refined. \n
 * @param[out] WORK
          WORK is REAL array, dimension (2*N) \n
          Workspace. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (2*N) \n
          Workspace. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot in the Sturm sequence for T. \n
 * @param[in] SPDIAM
          SPDIAM is REAL \n
          The spectral diameter of T. \n
 * @param[out] INFO
          INFO is INTEGER \n
          Error flag. \n

 *  * */
    template <typename T>
    void larrj(integer *n, T *d, T *e2, integer *ifirst, integer *ilast, T *rtol, integer *offset,
               T *w, T *werr, T *work, integer *iwork, T *pivmin, T *spdiam, integer *info)
    {
        larrj(n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam, info);
    }
    /** @}*/ // end of larrj

    /** @defgroup larrk larrk
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRK computes one eigenvalue of a symmetric tridiagonal matrix T to suitable
 accuracy

 * @details
 * \b Purpose:
    \verbatim
    LARRK computes one eigenvalue of a symmetric tridiagonal
    matrix T to suitable accuracy. This is an auxiliary code to be
    called from SSTEMR.

    To avoid overflow, the matrix must be scaled so that its
    largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and
 for greatest accuracy, it should not be much smaller than that.

    See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
    Matrix", Report CS41, Computer Science Dept., Stanford
    University, July 21, 1966.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the tridiagonal matrix T.  N >= 0. \n
 * @param[in] IW
          IW is INTEGER \n
          The index of the eigenvalues to be returned. \n
 * @param[in] GL
          GL is REAL \n
 * @param[in] GU
          GU is REAL \n
          An upper and a lower bound on the eigenvalue. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the tridiagonal matrix T. \n
 * @param[in] E2
          E2 is REAL array, dimension (N-1) \n
          The (n-1) squared off-diagonal elements of the tridiagonal matrix T. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot allowed in the Sturm sequence for T. \n
 * @param[in] RELTOL
          RELTOL is REAL \n
          The minimum relative width of an interval.  When an interval
          is narrower than RELTOL times the larger (in
          magnitude) endpoint, then it is considered to be
          sufficiently small, i.e., converged.  Note: this should
          always be at least radix*machine epsilon. \n
 * @param[out] W
          W is REAL \n
 * @param[out] WERR
          WERR is REAL \n
          The error bound on the corresponding eigenvalue approximation
          in W. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:       Eigenvalue converged \n
          = -1:      Eigenvalue did NOT converge  \n

 *  * */
    template <typename T>
    void larrk(integer *n, integer *iw, T *gl, T *gu, T *d, T *e2, T *pivmin, T *reltol, T *w,
               T *werr, integer *info)
    {
        larrk(n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info);
    }
    /** @}*/ // end of larrk

    /** @defgroup larrv larrv
     * @ingroup SYM
     * @{
     */
    /*! @brief LARRV computes the eigenvectors of the tridiagonal matrix T = L D LT    \n
     given L, D and the eigenvalues of L D LT
 * @details
 * \b Purpose:
    \verbatim
    LARRV computes the eigenvectors of the tridiagonal matrix
    T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
    The input eigenvalues should have been computed by LARRE.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 0. \n
 * @param[in] VL
          VL is REAL \n
          Lower bound of the interval that contains the desired
          eigenvalues. VL < VU. Needed to compute gaps on the left or right
          end of the extremal eigenvalues in the desired RANGE. \n
 * @param[in] VU
          VU is REAL \n
          Upper bound of the interval that contains the desired
          eigenvalues. VL < VU.  \n
          Note: VU is currently not used by this implementation of SLARRV, VU is
          passed to SLARRV because it could be used compute gaps on the right end
          of the extremal eigenvalues. However, with not much initial accuracy in
          LAMBDA and VU, the formula can lead to an overestimation of the right gap
          and thus to inadequately early RQI 'convergence'. This is currently
          prevented this by forcing a small right gap. And so it turns out that VU
          is currently not used by this implementation of SLARRV. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the N diagonal elements of the diagonal matrix D.
          On exit, D may be overwritten. \n
 * @param[in,out] L
          L is REAL array, dimension (N) \n
          On entry, the (N-1) subdiagonal elements of the unit
          bidiagonal matrix L are in elements 1 to N-1 of L
          (if the matrix is not split.) At the end of each block
          is stored the corresponding shift as given by SLARRE. \n
          On exit, L is overwritten. \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot allowed in the Sturm sequence. \n
 * @param[in] ISPLIT
          ISPLIT is INTEGER array, dimension (N) \n
          The splitting points, at which T breaks up into blocks.
          The first block consists of rows/columns 1 to
          ISPLIT(1), the second of rows/columns ISPLIT(1)+1
          through ISPLIT(2), etc. \n
 * @param[in] M
          M is INTEGER \n
          The total number of input eigenvalues.  0 <= M <= N. \n
 * @param[in] DOL
          DOL is INTEGER \n
 * @param[in] DOU
          DOU is INTEGER \n
          If the user wants to compute only selected eigenvectors from all
          the eigenvalues supplied, he can specify an index range DOL:DOU.
          Or else the setting DOL=1, DOU=M should be applied.
          Note that DOL and DOU refer to the order in which the eigenvalues
          are stored in W. \n
          If the user wants to compute only selected eigenpairs, then
          the columns DOL-1 to DOU+1 of the eigenvector space Z contain the
          computed eigenvectors. All other columns of Z are set to zero. \n
 * @param[in] MINRGP
          MINRGP is REAL \n
 * @param[in] RTOL1
          RTOL1 is REAL \n
 * @param[in] RTOL2
          RTOL2 is REAL \n
          Parameters for bisection. \n
          An interval [LEFT,RIGHT] has converged if
          RIGHT-LEFT < MAX(RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|)) \n
 * @param[in,out] W
          W is REAL array, dimension (N) \n
          The first M elements of W contain the APPROXIMATE eigenvalues for
          which eigenvectors are to be computed.  The eigenvalues
          should be grouped by split-off block and ordered from
          smallest to largest within the block (The output array
          W from SLARRE is expected here). Furthermore, they are with
          respect to the shift of the corresponding root representation
          for their block. On exit, W holds the eigenvalues of the
          UNshifted matrix. \n
 * @param[in,out] WERR
          WERR is REAL array, dimension (N) \n
          The first M elements contain the semiwidth of the uncertainty
          interval of the corresponding eigenvalue in W \n
 * @param[in,out] WGAP
          WGAP is REAL array, dimension (N) \n
          The separation from the right neighbor eigenvalue in W. \n
 * @param[in] IBLOCK
          IBLOCK is INTEGER array, dimension (N) \n
          The indices of the blocks (submatrices) associated with the
          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue
          W(i) belongs to the first block from the top, =2 if W(i)
          belongs to the second block, etc. \n
 * @param[in] INDEXW
          INDEXW is INTEGER array, dimension (N) \n
          The indices of the eigenvalues within each block (submatrix);
          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the
          i-th eigenvalue W(i) is the 10-th eigenvalue in the second block. \n
 * @param[in] GERS
          GERS is REAL array, dimension (2*N) \n
          The N Gerschgorin intervals (the i-th Gerschgorin interval
          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should
          be computed from the original UNshifted matrix. \n
 * @param[out] Z
          Z is REAL array, dimension (LDZ, fla_max(1,M)) \n
          If INFO = 0, the first M columns of Z contain the
          orthonormal eigenvectors of the matrix T
          corresponding to the input eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          Note: the user must ensure that at least fla_max(1,M) columns are
          supplied in the array Z. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(1,N). \n
 * @param[out] ISUPPZ
          ISUPPZ is INTEGER array, dimension (2*max(1,M)) \n
          The support of the eigenvectors in Z, i.e., the indices
          indicating the nonzero elements in Z. The I-th eigenvector
          is nonzero only in elements ISUPPZ(2*I-1) through
          ISUPPZ(2*I). \n
 * @param[out] WORK
          WORK is REAL array, dimension (12*N) \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (7*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          > 0:  A problem occurred in SLARRV. \n
          < 0:  One of the called subroutines signaled an internal problem.
                Needs inspection of the corresponding parameter IINFO
                for further information. \n
          =-1:  Problem in SLARRB when refining a child's eigenvalues. \n
          =-2:  Problem in SLARRF when computing the RRR of a child.
                When a child is inside a tight cluster, it can be difficult
                to find an RRR. A partial remedy from the user's point of
                view is to make the parameter MINRGP smaller and recompile.
                However, as the orthogonality of the computed vectors is
                proportional to 1/MINRGP, the user should be aware that
                he might be trading in precision when he decreases MINRGP. \n
          =-3:  Problem in SLARRB when refining a single eigenvalue
                after the Rayleigh correction was rejected. \n
          = 5:  The Rayleigh Quotient Iteration failed to converge to
                full accuracy in MAXITR steps.  \n

 *  * */
    template <typename T>
    void larrv(integer *n, T *vl, T *vu, T *d, T *l, T *pivmin, integer *isplit, integer *m,
               integer *dol, integer *dou, T *minrgp, T *rtol1, T *rtol2, T *w, T *werr, T *wgap,
               integer *iblock, integer *indexw, T *gers, T *z, integer *ldz, integer *isuppz,
               T *work, integer *iwork, integer *info)
    {
        larrv(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap,
              iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
    }
    template <typename T, typename Ta>
    void larrv(integer *n, Ta *vl, Ta *vu, Ta *d, Ta *l, Ta *pivmin, integer *isplit, integer *m,
               integer *dol, integer *dou, Ta *minrgp, Ta *rtol1, Ta *rtol2, Ta *w, Ta *werr,
               Ta *wgap, integer *iblock, integer *indexw, Ta *gers, T *z, integer *ldz,
               integer *isuppz, Ta *work, integer *iwork, integer *info)
    {
        larrv(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap,
              iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
    }
    /** @}*/ // end of larrv

    /** @defgroup lar1v lar1v
     * @ingroup SYM
     * @{
     */
    /*! @brief LAR1V computes the (scaled) r-th column of the inverse of the \n
     submatrix in rows b1 through bn of the tridiagonal matrix LDLT - Î»I
 * @details
 * \b Purpose:
    \verbatim
    LAR1V computes the (scaled) r-th column of the inverse of
    the sumbmatrix in rows B1 through BN of the tridiagonal matrix
    L D L**T - sigma I. When sigma is close to an eigenvalue, the
    computed vector is an accurate eigenvector. Usually, r corresponds
    to the index where the eigenvector is largest in magnitude.
    The following steps accomplish this computation :
    (a) Stationary qd transform, L D L**T - sigma I = L(+) D(+) L(+)**T,
    (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
    (c) Computation of the diagonal elements of the inverse of
        L D L**T - sigma I by combining the above transforms, and choosing
        r as the index where the diagonal of the inverse is (one of the)
        largest in magnitude.
    (d) Computation of the (scaled) r-th column of the inverse using the
        twisted factorization obtained by combining the top part of the
        the stationary and the bottom part of the progressive transform.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix L D L**T. \n
 * @param[in] B1
          B1 is INTEGER \n
          First index of the submatrix of L D L**T. \n
 * @param[in] BN
          BN is INTEGER \n
          Last index of the submatrix of L D L**T. \n
 * @param[in] LAMBDA
          LAMBDA is REAL \n
          The shift. In order to compute an accurate eigenvector,
          LAMBDA should be a good approximation to an eigenvalue
          of L D L**T. \n
 * @param[in] L
          L is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of the unit bidiagonal matrix
          L, in elements 1 to N-1. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the diagonal matrix D. \n
 * @param[in] LD
          LD is REAL array, dimension (N-1) \n
          The n-1 elements L(i)*D(i). \n
 * @param[in] LLD
          LLD is REAL array, dimension (N-1) \n
          The n-1 elements L(i)*L(i)*D(i). \n
 * @param[in] PIVMIN
          PIVMIN is REAL \n
          The minimum pivot in the Sturm sequence. \n
 * @param[in] GAPTOL
          GAPTOL is REAL \n
          Tolerance that indicates when eigenvector entries are negligible
          w.r.t. their contribution to the residual. \n
 * @param[in,out] Z
          Z is REAL array, dimension (N) \n
          On input, all entries of Z must be set to 0.
          On output, Z contains the (scaled) r-th column of the
          inverse. The scaling is such that Z(R) equals 1. \n
 * @param[in] WANTNC
          WANTNC is LOGICAL \n
          Specifies whether NEGCNT has to be computed. \n
 * @param[out] NEGCNT
          NEGCNT is INTEGER \n
          If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin
          in the  matrix factorization L D L**T, and NEGCNT = -1 otherwise. \n
 * @param[out] ZTZ
          ZTZ is REAL \n
          The square of the 2-norm of Z. \n
 * @param[out] MINGMA
          MINGMA is REAL \n
          The reciprocal of the largest (in magnitude) diagonal
          element of the inverse of L D L**T - sigma I. \n
 * @param[in,out] R
          R is INTEGER \n
          The twist index for the twisted factorization used to
          compute Z. \n
          On input, 0 <= R <= N. If R is input as 0, R is set to
          the index where (L D L**T - sigma I)^{-1} is largest
          in magnitude. If 1 <= R <= N, R is unchanged.
          On output, R contains the twist index used to compute Z.
          Ideally, R designates the position of the maximum entry in the
          eigenvector. \n
 * @param[out] ISUPPZ
          ISUPPZ is INTEGER array, dimension (2) \n
          The support of the vector in Z, i.e., the vector Z is
          nonzero only in elements ISUPPZ(1) through ISUPPZ(2). \n
 * @param[out] NRMINV
          NRMINV is REAL \n
          NRMINV = 1/SQRT(ZTZ) \n
 * @param[out] RESID
          RESID is REAL \n
          The residual of the FP vector.
          RESID = ABS(MINGMA)/SQRT(ZTZ) \n
 * @param[out] RQCORR
          RQCORR is REAL \n
          The Rayleigh Quotient correction to LAMBDA.
          RQCORR = MINGMA*TMP \n
 * @param[out] WORK
          WORK is REAL array, dimension (4*N)  \n

 *  * */
    template <typename T>
    void lar1v(integer *n, integer *b1, integer *bn, T *lambda, T *d, T *l, T *ld, T *lld,
               T *pivmin, T *gaptol, T *z, logical *wantnc, integer *negcnt, T *ztz, T *mingma,
               integer *r, integer *isuppz, T *nrminv, T *resid, T *rqcorr, T *work)
    {
        lar1v(n, b1, bn, lambda, d, l, d, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r,
              isuppz, nrminv, resid, rqcorr, work);
    }
    template <typename T, typename Ta>
    void lar1v(integer *n, integer *b1, integer *bn, Ta *lambda, Ta *d, Ta *l, Ta *ld, Ta *lld,
               Ta *pivmin, Ta *gaptol, T *z, logical *wantnc, integer *negcnt, Ta *ztz, Ta *mingma,
               integer *r, integer *isuppz, Ta *nrminv, Ta *resid, Ta *rqcorr, Ta *work)
    {
        lar1v(n, b1, bn, lambda, d, l, d, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r,
              isuppz, nrminv, resid, rqcorr, work);
    }
    /** @}*/ // end of lar1v
    /** @}*/ // end of SYM

    /** @defgroup HER Hermitian eigenvalues
     * @ingroup Eig
     * @{
     */

    /** @defgroup heev heev
     * @ingroup HER
     * @{
     */
    /*! @brief HEEV computes the eigenvalues and, optionally eigenvectors for complex HE matrices
               using eigenvalue decomposition (QR algorithm)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Eigenvalue decomposition (QR algorithm).
        Computation of all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix a.
    \endverbatim

    * @param[in] JOBS
              JOBZ is CHARACTER*1 \n
              = 'N':  Compute eigenvalues only; \n
              = 'V':  Compute eigenvalues and eigenvectors. \n
    * @param[in] UPLO
              UPLO is CHARACTER*1 \n
              = 'U':  Upper triangle of A is stored; \n
              = 'L':  Lower triangle of A is stored. \n
    * @param[in] N
              N is INTEGER \n
              The order of the matrix A.  N >= 0. \n
    * @param[in,out] A
              A is COMPLEX array, dimension (LDA, N) \n
              On entry, the Hermitian matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of A contains the
              upper triangular part of the matrix A.  If UPLO = 'L',
              the leading N-by-N lower triangular part of A contains
              the lower triangular part of the matrix A. \n
              On exit, if JOBZ = 'V', then if INFO = 0, A contains the
              orthonormal eigenvectors of the matrix A. \n
              If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
              or the upper triangle (if UPLO='U') of A, including the
              diagonal, is destroyed. \n
    * @param[in] LDA
              LDA is INTEGER \n
              The leading dimension of the array A.  LDA >= max(1,N). \n
    * @param[out] W
              w is REAL array, dimension (N) \n
              If INFO = 0, the eigenvalues in ascending order. \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The length of the array WORK.  LWORK >= max(1,2*N-1). \n
              For optimal efficiency, LWORK >= (NB+1)*N,
              where NB is the blocksize for CHETRD returned by ILAENV. \n

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out] RWORK
          RWORK is REAL array, dimension (max(1, 3*N-2)) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  if INFO = i, the algorithm failed to converge; i
                    off-diagonal elements of an intermediate tridiagonal
                    form did not converge to zero. \n

    *     *  */
    template< typename T, typename Ta >
    void heev(char* jobz, char* uplo, integer* n, T* a, integer* lda, Ta*  w, T* work, integer* lwork,
              Ta* rwork, integer* info)
    {
         heev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
    }
    /** @}*/ // end of heev

    /** @defgroup heevd heevd
     * @ingroup HER
     * @{
     */
    /*! @brief HEEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices
               using divide and conquer algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of the eigenvalues and, optionally, the left and/or right eigenvectors for HE
        matrices. HEEVD computes all eigenvalues and, optionally, eigenvectors of a comple Hermitian
        matrix A. If eigenvectors are desired, it uses a divide and conquer algorithm.

    \endverbatim

    * @param[in] JOBZ
              JOBZ is CHARACTER*1 \n
              = 'N':  Compute eigenvalues only; \n
              = 'V':  Compute eigenvalues and eigenvectors. \n
    * @param[in] UPLO
              UPLO is CHARACTER*1 \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] N
              N is INTEGER \n
              The order of the matrix A.  N >= 0. \n
    * @param[in,out] A
              A is COMPLEX array, dimension (LDA, N) \n
              On entry, the Hermitian matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of A contains the
              upper triangular part of the matrix A.  If UPLO = 'L',
              the leading N-by-N lower triangular part of A contains
              the lower triangular part of the matrix A. \n
              On exit, if JOBZ = 'V', then if INFO = 0, A contains the
              orthonormal eigenvectors of the matrix A. \n
              If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
              or the upper triangle (if UPLO='U') of A, including the
              diagonal, is destroyed. \n
    * @param[in] LDA
              LDA is INTEGER \n
              The leading dimension of the array A.  LDA >= max(1,N). \n
    * @param[out] W
              W is COMPLEX array, dimension (N) \n
              If INFO = 0, the eigenvalues in ascending order. \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              If N <= 1,               LWORK must be at least 1. \n
              If JOBZ = 'N' and N > 1, LWORK must be at least N+1. \n
              If JOBZ = 'V' and N > 1, LWORK must be at least
                                                    2*N + N**2. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal sizes of the WORK, RWORK and IWORK
              arrays, returns these values as the first entries of the WORK, 
              RWORK and IWORK arrays, and no error message related to LWORK or
              LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	RWORK
              RWORK is REAL array, dimension (MAX(1,LRWORK)) \n
              On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. \n
    * @param[in]	LRWORK
              LRWORK is INTEGER \n
              The dimension of the array RWORK. \n
              If N <= 1,                LRWORK must be at least 1. \n
              If JOBZ  = 'N' and N > 1, LRWORK must be at least N. \n
              If JOBZ  = 'V' and N > 1, LRWORK must be at least 1 + 5*N + 2*N**2. \n
 \n
              If LRWORK = -1, then a workspace query is assumed; the
              routine only calculates the optimal sizes of the WORK, RWORK and
              IWORK arrays, returns these values as the first entries of
              the WORK, RWORK and IWORK arrays, and no error message related to
              LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	IWORK
              IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
              On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. \n
    * @param[in]	LIWORK
              LIWORK is INTEGER \n
              The dimension of the array IWORK. \n
              If N <= 1,                LIWORK must be at least 1. \n
              If JOBZ  = 'N' and N > 1, LIWORK must be at least 1. \n
              If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. \n
 \n
              If LIWORK = -1, then a workspace query is assumed; the
              routine only calculates the optimal sizes of the WORK, RWORK and
              IWORK arrays, returns these values as the first entries of
              the WORK, RWORK and IWORK arrays, and no error message related to
              LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
                    to converge; i off-diagonal elements of an intermediate
                    tridiagonal form did not converge to zero;
                    if INFO = i and JOBZ = 'V', then the algorithm failed
                    to compute an eigenvalue while working on the submatrix
                    lying in rows and columns INFO/(N+1) through
                    mod(INFO,N+1). \n

    *     *  */
   template< typename T, typename Ta >
   void heevd(char* jobz, char* uplo, integer* n, T* a, integer* lda, Ta* w, T* work, 
              integer* lwork, Ta* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
   {
      heevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
   }
    /** @}*/ // end of heevd

    /** @defgroup heevr heevr
     * @ingroup HER
     * @{
     */
    /*! @brief HEEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices
               using Hermitian eigenvalue decomposition (MRRR)
    *

    * @details
    * \b Purpose:
    * \verbatim
         Hermitian eigenvalue decomposition (MRRR).
         HEEVR computes selected eigenvalues and, optionally, eigenvectors
         of a complex Hermitian matrix A. Eigenvalues and eigenvectors can be
         selected by specifying either a range of values or a range of indices
         for the desired eigenvalues. Invocations with different choices for
         these parameters may result in the computation of slightly different
         eigenvalues and/or eigenvectors for the same matrix. The reason for
         this behavior is that there exists a variety of algorithms (each
         performing best for a particular set of options) with HEEVR
         attempting to select the best based on the various parameters. In all
         cases, the computed values are accurate within the limits of finite
         precision arithmetic.

         HEEVR first reduces the matrix A to tridiagonal form T with a call
         to CHETRD.  Then, whenever possible, HEEVR calls CSTEMR to compute
         the eigenspectrum using Relatively Robust Representations.  CSTEMR
         computes eigenvalues by the dqds algorithm, while orthogonal
         eigenvectors are computed from various  L D L^T representations
         (also known as Relatively Robust Representations). Gram-Schmidt
         orthogonalization is avoided as far as possible. More specifically,
         the various steps of the algorithm are as follows.

         For each unreduced block (submatrix) of T,
         (a) Compute T - sigma I  = L D L^T, so that L and D
         define all the wanted eigenvalues to high relative accuracy.
         This means that small relative changes in the entries of D and L
         cause only small relative changes in the eigenvalues and
         eigenvectors. The standard (unfactored) representation of the
         tridiagonal matrix T does not have this property in general.
         (b) Compute the eigenvalues to suitable accuracy.
         If the eigenvectors are desired, the algorithm attains full
         accuracy of the computed eigenvalues only right before
         the corresponding vectors have to be computed, see steps c) and d).
         (c) For each cluster of close eigenvalues, select a new
         shift close to the cluster, find a new factorization, and refine
         the shifted eigenvalues to suitable accuracy.
         (d) For each eigenvalue with a large enough relative separation compute
         the corresponding eigenvector by forming a rank revealing twisted
         factorization. Go back to (c) for any clusters that remain.

         The desired accuracy of the output can be specified by the input
         parameter ABSTOL.

         For more details, see CSTEMR's documentation and:
         - Inderjit S. Dhillon and Beresford N. Parlett:
         Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
         - Inderjit Dhillon and Beresford Parlett:  SIAM Journal on Matrix Analysis and Applications, Vol. 25,
         2004.  Also LAPACK Working Note 154.
         - Inderjit Dhillon: ,
         Computer Science Division Technical Report No. UCB/CSD-97-971,
         UC Berkeley, May 1997.
    \endverbatim

    * @param[in] JOBZ
              JOBZ is CHARACTER*1 \n
              = 'N':  Compute eigenvalues only; \n
              = 'V':  Compute eigenvalues and eigenvectors. \n
    * @param[in] RANGE
              RANGE is CHARACTER*1 \n
              = 'A': all eigenvalues will be found. \n
              = 'V': all eigenvalues in the half-open interval (VL,VU]
              will be found. \n
              = 'I': the IL-th through IU-th eigenvalues will be found. \n
              For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and
              CSTEIN are called \n
    * @param[in] UPLO
              UPLO is CHARACTER*1 \n
              = 'U':  Upper triangle of A is stored; \n
              = 'L':  Lower triangle of A is stored. \n
    * @param[in] N
              N is INTEGER \n
              The order of the matrix A.  N >= 0. \n
    * @param[in,out] A
              A is COMPLEX array, dimension (LDA, N) \n
              On entry, the Hermitian matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of a contains the
              upper triangular part of the matrix A.  If UPLO = 'L',
              the leading N-by-N lower triangular part of A contains
              the lower triangular part of the matrix A. \n
              On exit, the lower triangle (if UPLO='L') or the upper
              triangle (if UPLO='U') of A, including the diagonal, is
              destroyed. \n
    * @param[in] LDA
              LDA is integer* \n
              The leading dimension of the array A.  LDA >= MAX(1,N). \n
    * @param[in] VL
              VL is REAL \n
              If RANGE='V', the lower bound of the interval to
              be searched for eigenvalues. VL < VU. \n
              Not referenced if range = 'A' or 'I'. \n
    * @param[in] VU
				  VU is REAL \n
				  If RANGE='V', the upper bound of the interval to
				  be searched for eigenvalues. VL < VU. \n
				  Not referenced if RANGE = 'A' or 'I'. \n
    * @param[in] IL
				  IL is INTEGER \n
				  If RANGE='I', the index of the
				  smallest eigenvalue to be returned. \n
				  1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
				  Not referenced if RANGE = 'A' or 'V'. \n
    * @param[in] IU
				  IU is INTEGER \n
				  If RANGE='I', the index of the
				  largest eigenvalue to be returned. \n
				  1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. \n
				  Not referenced if RANGE = 'A' or 'V'. \n
    * @param[in] ABSTOL
				  ABSTOL is REAL \n
				  The absolute error tolerance for the eigenvalues. \n
				  An approximate eigenvalue is accepted as converged
				  when it is determined to lie in an interval [a,b]
				  of width less than or equal to \n
				  
				  ABSTOL + EPS *   max( |a|,|b| ) , \n
				  
				  where EPS is the machine precision.  If ABSTOL is less than
				  or equal to zero, then  EPS*|T|  will be used in its place,
				  where |T| is the 1-norm of the tridiagonal matrix obtained
				  by reducing A to tridiagonal form. \n
				  
				  See  by Demmel and
				  Kahan, LAPACK Working Note #3. \n
				  
				  If high relative accuracy is important, set ABSTOL to
				  SLAMCH( 'Safe minimum' ).  Doing so will guarantee that
				  eigenvalues are computed to high relative accuracy when
				  possible in future releases.  The current code does not
				  make any guarantees about high relative accuracy, but
				  future releases will. See J. Barlow and J. Demmel,
				  , LAPACK Working Note #7, for a discussion
				  of which matrices define their eigenvalues to high relative
				  accuracy. \n
    * @param[out] M
				  M is INTEGER \n
				  The total number of eigenvalues found.  0 <= M <= N. \n
				  If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. \n
    * @param[out] W
				  W is REAL array, dimension (N) \n
				  The first M elements contain the selected eigenvalues in
				  ascending order. \n
    * @param[out] Z
				  Z is COMPLEX array, dimension (LDZ, max(1,M)) \n
				  If JOBZ = 'V', then if INFO = 0, the first M columns of Z
				  contain the orthonormal eigenvectors of the matrix A
				  corresponding to the selected eigenvalues, with the i-th
				  column of Z holding the eigenvector associated with W(i). \n
				  If JOBZ = 'N', then Z is not referenced. \n
				  Note: the user must ensure that at least max(1,M) columns are
				  supplied in the array Z; if RANGE = 'V', the exact value of M
				  is not known in advance and an upper bound must be used.
				  Supplying N columns is always safe. \n
    * @param[in] LDZ
				  LDZ is INTEGER \n
				  The leading dimension of the array Z.  LDZ >= 1, and if
				  JOBZ = 'V', LDZ >= max(1,N). \n
    * @param[out] ISUPPZ
				  ISUPPZ is INTEGER array, dimension ( 2*max(1,M) ) \n
				  The support of the eigenvectors in Z, i.e., the indices
				  indicating the nonzero elements in Z. The i-th eigenvector
				  is nonzero only in elements ISUPPZ( 2*i-1 ) through
				  ISUPPZ( 2*i ). This is an output of CSTEMR (tridiagonal
				  matrix). The support of the eigenvectors of A is typically
				  1:N because of the unitary transformations applied by CUNMTR. \n
				  Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 \n
    * @param[out]	WORK
				  WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
				  On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
				  LWORK is INTEGER \n
				  The length of the array WORK. \n
				  If N <= 1, LWORK >= 1, else LWORK >= 2*N. \n
				  For optimal efficiency, LWORK >= (NB+1)*N,
				  where NB is the max of the blocksize for CHETRD and for
				  CUNMTR as returned by ILAENV. \n \n
				  
				  If LWORK = -1, then a workspace query is assumed; the routine
				  only calculates the optimal sizes of the WORK, RWORK and
				  IWORK arrays, returns these values as the first entries of
				  the WORK, RWORK and IWORK arrays, and no error message
				  related to LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	RWORK
				  RWORK is REAL array, dimension (MAX(1,LRWORK)) \n
				  On exit, if INFO = 0, RWORK(1) returns the optimal
				  (and minimal) LRWORK. \n
    * @param[in]	LRWORK
				  LRWORK is INTEGER \n
				  The length of the array RWORK. \n
				  If N <= 1, LRWORK >= 1, else LRWORK >= 24*N. \n \n
				  
				  If LRWORK = -1, then a workspace query is assumed; the
				  routine only calculates the optimal sizes of the WORK, RWORK
				  and IWORK arrays, returns these values as the first entries
				  of the WORK, RWORK and IWORK arrays, and no error message
				  related to LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	IWORK
				  IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
				  On exit, if INFO = 0, IWORK(1) returns the optimal
				  (and minimal) LIWORK. \n
    * @param[in]	LIWORK
				  LIWORK is INTEGER \n
				  The dimension of the array IWORK. \n
				  If N <= 1, LIWORK >= 1, else LIWORK >= 10*N. \n \n
				  
				  If LIWORK = -1, then a workspace query is assumed; the
				  routine only calculates the optimal sizes of the WORK, RWORK
				  and IWORK arrays, returns these values as the first entries
				  of the WORK, RWORK and IWORK arrays, and no error message
				  related to LWORK or LRWORK or LIWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  Internal error \n

    *     *  */
	 template< typename T, typename Ta >
	 void heevr(char* jobz, char* range, char* uplo, integer* n, T* a, integer* lda,
               Ta*  vl, Ta*  vu, integer* il, integer* iu, Ta*  abstol, integer* m, Ta* w,
               T* z, integer* ldz, integer* isuppz, T* work, integer* lwork, Ta* rwork, integer* lrwork, 
               integer* iwork, integer* liwork, integer* info)
	 {
	   heevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, 
            work, lwork, rwork, lrwork, iwork, liwork, info);
	 }
    /** @}*/ // end of heevr

    /** @}*/ // end of HER
    /** @}*/ // end of Eig
}
#endif