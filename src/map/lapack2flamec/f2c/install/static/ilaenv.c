/* ./ilaenv.f -- translated by f2c (version 20190311). You must link the
 resulting object file with libf2c: on Microsoft Windows system, link with
 libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install
 libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of
 the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
/* > \brief \b ILAENV */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ILAENV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 ) */
/* .. Scalar Arguments .. */
/* CHARACTER*( * ) NAME, OPTS */
/* INTEGER ISPEC, N1, N2, N3, N4 */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ILAENV is called from the LAPACK routines to choose problem-dependent */
/* > parameters for the local environment. See ISPEC for a description of */
/* > the parameters. */
/* > */
/* > ILAENV returns an INTEGER */
/* > if ILAENV >= 0: ILAENV returns the value of the parameter specified by
 * ISPEC */
/* > if ILAENV < 0: if ILAENV = -k, the k-th argument had an illegal value. */
/* > */
/* > This version provides a set of parameters which should give good, */
/* > but not optimal, performance on many of the currently available */
/* > computers. Users are encouraged to modify this subroutine to set */
/* > the tuning parameters for their particular machine using the option */
/* > and problem size information in the arguments. */
/* > */
/* > This routine will not function correctly if it is converted to all */
/* > lower case. Converting it to all upper case is allowed. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ISPEC */
/* > \verbatim */
/* > ISPEC is INTEGER */
/* > Specifies the parameter to be returned as the value of */
/* > ILAENV. */
/* > = 1: the optimal blocksize;
if this value is 1, an unblocked */
/* > algorithm will give the best performance. */
/* > = 2: the minimum block size for which the block routine */
/* > should be used;
if the usable block size is less than */
/* > this value, an unblocked routine should be used. */
/* > = 3: the crossover point (in a block routine, for N less */
/* > than this value, an unblocked routine should be used) */
/* > = 4: the number of shifts, used in the nonsymmetric */
/* > eigenvalue routines (DEPRECATED) */
/* > = 5: the minimum column dimension for blocking to be used;
 */
/* > rectangular blocks must have dimension at least k by m, */
/* > where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
/* > = 6: the crossover point for the SVD (when reducing an m by n */
/* > matrix to bidiagonal form, if fla_max(m,n)/min(m,n) exceeds */
/* > this value, a QR factorization is used first to reduce */
/* > the matrix to a triangular form.) */
/* > = 7: the number of processors */
/* > = 8: the crossover point for the multishift QR method */
/* > for nonsymmetric eigenvalue problems (DEPRECATED) */
/* > = 9: maximum size of the subproblems at the bottom of the */
/* > computation tree in the divide-and-conquer algorithm */
/* > (used by xGELSD and xGESDD) */
/* > =10: ieee infinity and NaN arithmetic can be trusted not to trap */
/* > =11: infinity arithmetic can be trusted not to trap */
/* > 12 <= ISPEC <= 17: */
/* > xHSEQR or related subroutines, */
/* > see IPARMQ for detailed explanation */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* > NAME is CHARACTER*(*) */
/* > The name of the calling subroutine, in either upper case or */
/* > lower case. */
/* > \endverbatim */
/* > */
/* > \param[in] OPTS */
/* > \verbatim */
/* > OPTS is CHARACTER*(*) */
/* > The character options to the subroutine NAME, concatenated */
/* > into a single character string. For example, UPLO = 'U', */
/* > TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/* > be specified as OPTS = 'UTN'. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* > N1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* > N2 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N3 */
/* > \verbatim */
/* > N3 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N4 */
/* > \verbatim */
/* > N4 is INTEGER */
/* > Problem dimensions for the subroutine NAME;
these may not all */
/* > be required. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup ilaenv */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The following conventions have been used when calling ILAENV from the */
/* > LAPACK routines: */
/* > 1) OPTS is a concatenation of all of the character options to */
/* > subroutine NAME, in the same order that they appear in the */
/* > argument list for NAME, even if they are not used in determining */
/* > the value of the parameter specified by ISPEC. */
/* > 2) The problem dimensions N1, N2, N3, N4 are specified in the order */
/* > that they appear in the argument list for NAME. N1 is used */
/* > first, N2 second, and so on, and unused problem dimensions are */
/* > passed a value of -1. */
/* > 3) The parameter value returned by ILAENV is checked for validity in */
/* > the calling subroutine. For example, ILAENV is used to retrieve */
/* > the optimal blocksize for STRTRI as follows: */
/* > */
/* > NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/* > IF( NB.LE.1 ) NB = MAX( 1, N ) */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/** Generated wrapper function */
aocl_int_t ilaenv_(aocl_int_t *ispec, char *name__, char *opts, aocl_int_t *n1, aocl_int_t *n2, aocl_int_t *n3, aocl_int_t *n4)
{
#if FLA_ENABLE_ILP64
    return aocl_lapack_ilaenv(ispec, name__, opts, n1, n2, n3, n4);
#else
    aocl_int64_t ispec_64 = *ispec;
    aocl_int64_t n1_64 = *n1;
    aocl_int64_t n2_64 = *n2;
    aocl_int64_t n3_64 = *n3;
    aocl_int64_t n4_64 = *n4;

    return aocl_lapack_ilaenv(&ispec_64, name__, opts, &n1_64, &n2_64, &n3_64, &n4_64);
#endif
}
