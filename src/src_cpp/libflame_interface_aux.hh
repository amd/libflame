/******************************************************************************
 * Copyright (C) 2021-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface_aux.hh
 *  libflame_interfac_auxe.hh defines all the Auxilliary routines for Libflame CPP templated public
 *  interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_AUX_HH
#define LIBFLAME_INTERFACE_AUX_HH

#include "libflame.hh"

namespace libflame
{
     /** @defgroup Auxs Auxilliary routines
     * @ingroup AOCL_LAPACK
     * @{
     */

    /** @defgroup Parameters Parameters
     * @ingroup Auxs
     * @{
     */

    /** @defgroup lamch lamch
     * @ingroup Parameters
     * @{
     */
    /*! @brief LAMCH determines single precision machine parameters

 * @details
 * \b Purpose:
    \verbatim
     LAMCH determines single precision machine parameters.
    \endverbatim

  * @param[in] CMACH
           CMACH is CHARACTER*1 \n
           Specifies the value to be returned by SLAMCH: \n
           = 'E' or 'e',   SLAMCH := eps \n
           = 'S' or 's ,   SLAMCH := sfmin \n
           = 'B' or 'b',   SLAMCH := base \n
           = 'P' or 'p',   SLAMCH := eps*base \n
           = 'N' or 'n',   SLAMCH := t \n
           = 'R' or 'r',   SLAMCH := rnd \n
           = 'M' or 'm',   SLAMCH := emin \n
           = 'U' or 'u',   SLAMCH := rmin \n
           = 'L' or 'l',   SLAMCH := emax \n
           = 'O' or 'o',   SLAMCH := rmax \n
           where \n
           eps   = relative machine precision \n
           sfmin = safe minimum, such that 1/sfmin does not overflow \n
           base  = base of the machine \n
           prec  = eps*base \n
           t     = number of (base) digits in the mantissa \n
           rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise \n
           emin  = minimum exponent before (gradual) underflow \n
           rmin  = underflow threshold - base**(emin-1) \n
           emax  = largest exponent before overflow \n
           rmax  = overflow threshold  - (base**emax)*(1-eps) \n

* @return T Template based return value of the function.
 * */
    template <typename T>
    T lamch(char *cmach)
    {
        if(typeid(T) == typeid(float))
        {
            return slamch_(cmach);
        }
        else
        {
            return dlamch_(cmach);
        }
    }
    /**@} */ // end of lamch

    /** @defgroup labad labad
     * @ingroup Parameters
     * @{
     */
    /*! @brief LABAD takes as input the values computed by LAMCH for underflow and  \n
     overflow, and   returns the square root of each of these values

 * @details
 * \b Purpose:
    \verbatim
    LABAD takes as input the values computed by SLAMCH for underflow and
    overflow, and   returns the square root of each of these values if the
    log of LARGE is sufficiently large.  This subroutine is intended to
    identify machines with a large exponent range, such as the Crays, and
    redefine the underflow and overflow limits to be the square roots of
    the values computed by SLAMCH.  This subroutine is needed because
    SLAMCH does not compensate for poor arithmetic in the upper half of
    the exponent range, as is found on a Cray.
    \endverbatim

 * @param[in,out] SMALL_VAL
          SMALL_VAL is REAL \n
          On entry, the underflow threshold as computed by SLAMCH.
          On exit, if LOG10(LARGE) is sufficiently large, the square
          root of SMALL, otherwise unchanged. \n
 * @param[in,out] LARGE
          LARGE is REAL \n
          On entry, the overflow threshold as computed by SLAMCH.
          On exit, if LOG10(LARGE) is sufficiently large, the square
          root of LARGE, otherwise unchanged. \n

 * */
    template <typename T>
    void labad(T *small_val, T *large)
    {
        labad(small_val, large);
    }
    /**@} */ // end of labad

    /** @defgroup ilaver ilaver
     * @ingroup Parameters
     * @{
     */
    /*! @brief ILAVER   returns the LAPACK version

 * @details
 * \b Purpose:
    \verbatim
    This subroutine   returns the LAPACK version.
    \endverbatim

 * @param[out] VERS_MAJOR
          VERS_MAJOR is INTEGER \n
          return the lapack major version \n
 * @param[out] VERS_MINOR
          VERS_MINOR is INTEGER \n
          return the lapack minor version from the major version \n
 * @param[out] VERS_PATCH
          VERS_PATCH is INTEGER \n
          return the lapack patch version from the minor version \n

 * */
    inline void ilaver(integer *vers_major, integer *vers_minor, integer *vers_patch__)
    {
        ilaver_(vers_major, vers_minor, vers_patch__);
    }
    /**@} */ // end of ilaver

    /** @defgroup iparam2stage iparam2stage
     * @ingroup Parameters
     * @{
     */
    /*! @brief IPARAM2STAGE program sets problem and machine dependent parameters

 * @details
 * \b Purpose:
    \verbatim
     This program sets problem and machine dependent parameters
     useful for xHETRD_2STAGE, xHETRD_HE2HB, xHETRD_HB2ST,
     xGEBRD_2STAGE, xGEBRD_GE2GB, xGEBRD_GB2BD
     and related subroutines for eigenvalue problems.
     It is called whenever ILAENV is called with 17 <= ISPEC <= 21.
     It is called whenever ILAENV2STAGE is called with 1 <= ISPEC <= 5
     with a direct conversion ISPEC + 16.
    \endverbatim

  * @param[in] ISPEC
           ISPEC is integer scalar \n
           ISPEC specifies which tunable parameter IPARAM2STAGE should
           return. \n
 \n
           ISPEC=17: the optimal blocksize nb for the reduction to
                     BAND
 \n
           ISPEC=18: the optimal blocksize ib for the eigenvectors
                     singular vectors update routine
 \n
           ISPEC=19: The length of the array that store the Housholder
                     representation for the second stage
                     Band to Tridiagonal or Bidiagonal
 \n
           ISPEC=20: The workspace needed for the routine in input.
 \n
           ISPEC=21: For future release. \n
  * @param[in] NAME
           NAME is character string \n
           Name of the calling subroutine \n
  * @param[in] OPTS
           OPTS is CHARACTER*(*) \n
           The character options to the subroutine NAME, concatenated
           into a single character string.  For example, UPLO = 'U',
           TRANS = 'T', and DIAG = 'N' for a triangular routine would
           be specified as OPTS = 'UTN'. \n
  * @param[in] NI
           NI is INTEGER which is the size of the matrix \n
  * @param[in] NBI
           NBI is INTEGER which is the used in the reduciton,  \n
           (e.g., the size of the band), needed to compute workspace
           and LHOUS2. \n
  * @param[in] IBI
           IBI is INTEGER which represent the IB of the reduciton,
           needed to compute workspace and LHOUS2. \n
  * @param[in] NXI
           NXI is INTEGER needed in the future release. \n

 * @return INTEGER Return value of the function.
 *  * */
    inline integer iparam2stage(integer *ispec, char *name, char *opts, integer *ni, integer *nbi,
                                integer *ibi, integer *nxi)
    {
        return iparam2stage_(ispec, name, opts, ni, nbi, ibi, nxi);
    }
    /**@} */ // end of iparam2stage

    /** @defgroup iladiag iladiag
     * @ingroup Parameters
     * @{
     */
    /*! @brief ILADIAG translated from a character string specifying if a matrix   \n
     has unit diagonal or not to the relevant BLAST-specified integer constant
 * @details
 * \b Purpose:
    \verbatim
    This subroutine translated from a character string specifying if a
    matrix has unit diagonal or not to the relevant BLAST-specified
    integer constant.

    ILADIAG   returns an INTEGER. If ILADIAG < 0, then the input is not a
    character indicating a unit or non-unit diagonal. Otherwise ILADIAG
    returns the constant value corresponding to DIAG.
    \endverbatim
 * @param[in] diag
          DIAG is character pointer. \n
          = N means non-unit diagonal. \n
          = U means unit diagonal. \n

 * @return INTEGER Return value of the function.
 *  * */
    inline integer iladiag(char *diag)
    {
        return iladiag_(diag);
    }
    /**@} */ // end of iladiag

    /** @defgroup ilatrans ilatrans
     * @ingroup Parameters
     * @{
     */
    /*! @brief ILATRANS translates from a character string specifying a transposition \n
     operation to the relevant BLAST-specified integer constant
 * @details
 * \b Purpose:
    \verbatim
    This subroutine translates from a character string specifying a
    transposition operation to the relevant BLAST-specified integer
    constant.

    ILATRANS   returns an INTEGER.  If ILATRANS < 0, then the input is not
    a character indicating a transposition operator.  Otherwise ILATRANS
      returns the constant value corresponding to TRANS
    \endverbatim
 * @param[in] trans
          trans is charater pointer. \n
          = N means no transposition \n
          = T means transposition \n
          = C means conjugate transposition \n

 * @return INTEGER Return value of the function.
 *  * */
    inline integer ilatrans(char *trans)
    {
        return ilatrans_(trans);
    }
    /**@} */ // end of ilatrans

    /** @defgroup ilaprec ilaprec
     * @ingroup Parameters
     * @{
     */
    /*! @brief ILAPREC translated from a character string specifying an intermediate \n
     precision to the relevant BLAST-specified integer constant.
 * @details
 * \b Purpose:
    \verbatim
    This subroutine translated from a character string specifying an
    intermediate precision to the relevant BLAST-specified integer
    constant.

    ILAPREC   returns an INTEGER.  If ILAPREC < 0, then the input is not a
    character indicating a supported intermediate precision.  Otherwise
    ILAPREC   returns the constant value corresponding to PREC.
    \endverbatim
 * @param[in] prec
          PREC is character pointer. \n
          = S means single precision \n
          = D means double precision \n
          = I means indigenous precision \n
          = X or E means extra precision. \n
          Other values returns invalid result. \n

 * @return INTEGER Return value of the function.
 *  * */
    inline integer ilaprec(char *prec)
    {
        return ilaprec_(prec);
    }
    /**@} */ // end of ilaprec

    /** @defgroup ilaenv2stage ilaenv2stage
     * @ingroup Parameters
     * @{
     */
    /*! @brief ILAENV2STAGE is called from the LAPACK routines to choose problem-dependent \n
     parameters for the local environment
 * @details
 * \b Purpose:
    \verbatim
    ILAENV2STAGE is called from the LAPACK routines to choose problem-dependent
    parameters for the local environment.  See ISPEC for a description of
    the parameters.
    It sets problem and machine dependent parameters useful for *_2STAGE and
    related subroutines.

    ILAENV2STAGE   returns an INTEGER
    if ILAENV2STAGE >= 0: ILAENV2STAGE   returns the value of the parameter
                          specified by ISPEC
    if ILAENV2STAGE < 0:  if ILAENV2STAGE = -k, the k-th argument had an
                          illegal value.

    This version provides a set of parameters which should give good,
    but not optimal, performance on many of the currently available
    computers for the 2-stage solvers. Users are encouraged to modify this
    subroutine to set the tuning parameters for their particular machine using
    the option and problem size information in the arguments.

    This routine will not function correctly if it is converted to all
    lower case.  Converting it to all upper case is allowed.
    \endverbatim

 * @param[in] ISPEC
          ISPEC is INTEGER \n
          Specifies the parameter to be   returned as the value of
          ILAENV2STAGE. \n
          = 1: the optimal blocksize nb for the reduction to BAND
 \n
          = 2: the optimal blocksize ib for the eigenvectors
               singular vectors update routine
 \n
          = 3: The length of the array that store the Housholder
               representation for the second stage
               Band to Tridiagonal or Bidiagonal
 \n
          = 4: The workspace needed for the routine in input.
 \n
          = 5: For future release. \n
 * @param[in] NAME
          NAME is CHARACTER*(*) \n
          The name of the calling subroutine, in either upper case or
          lower case. \n
 * @param[in] OPTS
          OPTS is CHARACTER*(*) \n
          The character options to the subroutine NAME, concatenated
          into a single character string.  For example, UPLO = 'U',
          TRANS = 'T', and DIAG = 'N' for a triangular routine would
          be specified as OPTS = 'UTN'. \n
 * @param[in] N1
          N1 is INTEGER \n
 * @param[in] N2
          N2 is INTEGER \n
 * @param[in] N3
          N3 is INTEGER \n
 * @param[in] N4
          N4 is INTEGER \n
          Problem dimensions for the subroutine NAME; these may not all
          be required.  \n

 * @return INTEGER Return value of the function.
 *  * */
    inline integer ilaenv2stage(integer *ispec, char *name, char *opts, integer *n1, integer *n2,
                                integer *n3, integer *n4)
    {
        return ilaenv2stage_(ispec, name, opts, n1, n2, n3, n4);
    }
    /**@} */ // end of ilaenv2stage

    /** @defgroup ilauplo ilauplo
     * @ingroup Parameters
     * @{
     */
    /*! @brief ILAUPLO translated from a character string specifying a upper- or \n
     lower-triangular matrix to the relevant BLAST-specified integer constant
 * @details
 * \b Purpose:
    \verbatim
    This subroutine translated from a character string specifying a
    upper- or lower-triangular matrix to the relevant BLAST-specified
    integer constant.

    ILAUPLO   returns an INTEGER.  If ILAUPLO < 0, then the input is not
    a character indicating an upper- or lower-triangular matrix.
    Otherwise ILAUPLO   returns the constant value corresponding to UPLO
    \endverbatim
 * @param[in] uplo
          uplo is charater pointer. \n
          = 'L' means lower triangular matrix \n
          = 'U' means upper triangular matrix \n
          Other values, returns invalid result. \n

 * @return INTEGER Return value of the function.
 *  * */
    inline integer ilauplo(char *uplo)
    {
        return ilauplo_(uplo);
    }
    /**@} */ // end of ilauplo

    /** @defgroup latrans la_trans
     * @ingroup Parameters
     * @{
     */
    /*! @brief CHLA_TRANSTYPE subroutine translates from a BLAST-specified integer constant \n
     to the character string specifying a transposition operation

 * @details
 * \b Purpose:
    \verbatim
    This subroutine translates from a BLAST-specified integer constant to
    the character string specifying a transposition operation.

    CHLA_TRANSTYPE   returns an CHARACTER*1.  If CHLA_TRANSTYPE is 'X',
    then input is not an integer indicating a transposition operator.
    Otherwise CHLA_TRANSTYPE   returns the constant value corresponding to
    TRANS.
    \endverbatim
 * */
    inline void hla_transtype(char *ret_val, integer *trans)
    {
        chla_transtype_(ret_val, trans);
    }
    /**@} */ // end of latrans
    /**@} */ // end of Parameters
    /**@} */ // end of Auxs
}
#endif