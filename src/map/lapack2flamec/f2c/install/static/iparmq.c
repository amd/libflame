/* ../netlib/v3.9.0/iparmq.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b IPARMQ */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download IPARMQ + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK ) */
/* .. Scalar Arguments .. */
/* INTEGER IHI, ILO, ISPEC, LWORK, N */
/* CHARACTER NAME*( * ), OPTS*( * ) */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This program sets problem and machine dependent parameters */
/* > useful for xHSEQR and related subroutines for eigenvalue */
/* > problems. It is called whenever */
/* > IPARMQ is called with 12 <= ISPEC <= 16 */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ISPEC */
/* > \verbatim */
/* > ISPEC is INTEGER */
/* > ISPEC specifies which tunable parameter IPARMQ should */
/* > return. */
/* > */
/* > ISPEC=12: (INMIN) Matrices of order nmin or less */
/* > are sent directly to xLAHQR, the implicit */
/* > double shift QR algorithm. NMIN must be */
/* > at least 11. */
/* > */
/* > ISPEC=13: (INWIN) Size of the deflation window. */
/* > This is best set greater than or equal to */
/* > the number of simultaneous shifts NS. */
/* > Larger matrices benefit from larger deflation */
/* > windows. */
/* > */
/* > ISPEC=14: (INIBL) Determines when to stop nibbling and */
/* > invest in an (expensive) multi-shift QR sweep. */
/* > If the aggressive early deflation subroutine */
/* > finds LD converged eigenvalues from an order */
/* > NW deflation window and LD > (NW*NIBBLE)/100, */
/* > then the next QR sweep is skipped and early */
/* > deflation is applied immediately to the */
/* > remaining active diagonal block. Setting */
/* > IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a */
/* > multi-shift QR sweep whenever early deflation */
/* > finds a converged eigenvalue. Setting */
/* > IPARMQ(ISPEC=14) greater than or equal to 100 */
/* > prevents TTQRE from skipping a multi-shift */
/* > QR sweep. */
/* > */
/* > ISPEC=15: (NSHFTS) The number of simultaneous shifts in */
/* > a multi-shift QR iteration. */
/* > */
/* > ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the */
/* > following meanings. */
/* > 0: During the multi-shift QR/QZ sweep, */
/* > blocked eigenvalue reordering, blocked */
/* > Hessenberg-triangular reduction, */
/* > reflections and/or rotations are not */
/* > accumulated when updating the */
/* > far-from-diagonal matrix entries. */
/* > 1: During the multi-shift QR/QZ sweep, */
/* > blocked eigenvalue reordering, blocked */
/* > Hessenberg-triangular reduction, */
/* > reflections and/or rotations are */
/* > accumulated, and matrix-matrix */
/* > multiplication is used to update the */
/* > far-from-diagonal matrix entries. */
/* > 2: During the multi-shift QR/QZ sweep, */
/* > blocked eigenvalue reordering, blocked */
/* > Hessenberg-triangular reduction, */
/* > reflections and/or rotations are */
/* > accumulated, and 2-by-2 block structure */
/* > is exploited during matrix-matrix */
/* > multiplies. */
/* > (If xTRMM is slower than xGEMM, then */
/* > IPARMQ(ISPEC=16)=1 may be more efficient than */
/* > IPARMQ(ISPEC=16)=2 despite the greater level of */
/* > arithmetic work implied by the latter choice.) */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* > NAME is CHARACTER string */
/* > Name of the calling subroutine */
/* > \endverbatim */
/* > */
/* > \param[in] OPTS */
/* > \verbatim */
/* > OPTS is CHARACTER string */
/* > This is a concatenation of the string arguments to */
/* > TTQRE. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > N is the order of the Hessenberg matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > It is assumed that H is already upper triangular */
/* > in rows and columns 1:ILO-1 and IHI+1:N. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The amount of workspace available. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup OTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Little is known about how best to choose these parameters. */
/* > It is possible to use different values of the parameters */
/* > for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR. */
/* > */
/* > It is probably best to choose different parameters for */
/* > different matrices and different parameters at different */
/* > times during the iteration, but this has not been */
/* > implemented --- yet. */
/* > */
/* > */
/* > The best choices of most of the parameters depend */
/* > in an ill-understood way on the relative execution */
/* > rate of xLAQR3 and xLAQR5 and on the nature of each */
/* > particular eigenvalue problem. Experiment may be the */
/* > only practical way to determine which choices are most */
/* > effective. */
/* > */
/* > Following is a list of default values supplied by IPARMQ. */
/* > These defaults may be adjusted in order to attain better */
/* > performance in any particular computational environment. */
/* > */
/* > IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point. */
/* > Default: 75. (Must be at least 11.) */
/* > */
/* > IPARMQ(ISPEC=13) Recommended deflation window size. */
/* > This depends on ILO, IHI and NS, the */
/* > number of simultaneous shifts returned */
/* > by IPARMQ(ISPEC=15). The default for */
/* > (IHI-ILO+1) <= 500 is NS. The default */
/* > for (IHI-ILO+1) > 500 is 3*NS/2. */
/* > */
/* > IPARMQ(ISPEC=14) Nibble crossover point. Default: 14. */
/* > */
/* > IPARMQ(ISPEC=15) Number of simultaneous shifts, NS. */
/* > a multi-shift QR iteration. */
/* > */
/* > If IHI-ILO+1 is ... */
/* > */
/* > greater than ...but less ... the */
/* > or equal to ... than default is */
/* > */
/* > 0 30 NS = 2+ */
/* > 30 60 NS = 4+ */
/* > 60 150 NS = 10 */
/* > 150 590 NS = ** */
/* > 590 3000 NS = 64 */
/* > 3000 6000 NS = 128 */
/* > 6000 infinity NS = 256 */
/* > */
/* > (+) By default matrices of this order are */
/* > passed to the implicit double shift routine */
/* > xLAHQR. See IPARMQ(ISPEC=12) above. These */
/* > values of NS are used only in case of a rare */
/* > xLAHQR failure. */
/* > */
/* > (**) The asterisks (**) indicate an ad-hoc */
/* > function increasing from 10 to 64. */
/* > */
/* > IPARMQ(ISPEC=16) Select structured matrix multiply. */
/* > (See ISPEC=16 above for details.) */
/* > Default: 3. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/** Generated wrapper function */
aocl_int_t iparmq_(aocl_int_t *ispec, char *name__, char *opts, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *lwork)
{
#if FLA_ENABLE_ILP64
    return aocl_lapack_iparmq(ispec, name__, opts, n, ilo, ihi, lwork);
#else
    aocl_int64_t ispec_64 = *ispec;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ilo_64 = *ilo;
    aocl_int64_t ihi_64 = *ihi;
    aocl_int64_t lwork_64 = *lwork;

    return aocl_lapack_iparmq(&ispec_64, name__, opts, &n_64, &ilo_64, &ihi_64, &lwork_64);
#endif
}

