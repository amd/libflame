/* ../netlib/ssbev.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b11 = 1.f;
static integer c__1 = 1;
/* > \brief <b> SSBEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for OTHER m atrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSBEV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbev.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbev.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbev.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KD, LDAB, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBEV computes all the eigenvalues and, optionally, eigenvectors of */
/* > a real symmetric band matrix A. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBZ */
/* > \verbatim */
/* > JOBZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only;
 */
/* > = 'V': Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is REAL array, dimension (LDAB, N) */
/* > On entry, the upper or lower triangle of the symmetric band */
/* > matrix A, stored in the first KD+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+kd). */
/* > */
/* > On exit, AB is overwritten by values generated during the */
/* > reduction to tridiagonal form. If UPLO = 'U', the first */
/* > superdiagonal and the diagonal of the tridiagonal matrix T */
/* > are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
/* > the diagonal and first subdiagonal of T are returned in the */
/* > first two rows of AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD + 1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* > eigenvectors of the matrix A, with the i-th column of Z */
/* > holding the eigenvector associated with W(i). */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1, and if */
/* > JOBZ = 'V', LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (fla_max(1,3*N-2)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the algorithm failed to converge;
i */
/* > off-diagonal elements of an intermediate tridiagonal */
/* > form did not converge to zero. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHEReigen */
/* ===================================================================== */
/* Subroutine */
void ssbev_(char *jobz, char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *w,
            real *z__, integer *ldz, real *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
             "ssbev inputs: jobz %c, uplo %c, n %" FLA_IS ", kd %" FLA_IS ", ldab %" FLA_IS
             ", ldz %" FLA_IS "",
             *jobz, *uplo, *n, *kd, *ldab, *ldz);
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real eps;
    integer inde;
    real anrm;
    integer imax;
    real rmin, rmax, sigma;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    logical lower, wantz;
    integer iscale;
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    extern real slansb_(char *, char *, integer *, integer *, real *, integer *, real *);
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *);
    integer indwrk;
    extern /* Subroutine */
        void
        ssbtrd_(char *, char *, integer *, integer *, real *, integer *, real *, real *, real *,
                integer *, real *, integer *),
        ssterf_(integer *, real *, real *, integer *);
    real smlnum;
    extern /* Subroutine */
        void
        ssteqr_(char *, integer *, real *, real *, real *, integer *, real *, integer *);
    /* -- LAPACK driver routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(lower || lsame_(uplo, "U", 1, 1)))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*kd < 0)
    {
        *info = -4;
    }
    else if(*ldab < *kd + 1)
    {
        *info = -6;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -9;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSBEV ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 1)
    {
        if(lower)
        {
            w[1] = ab[ab_dim1 + 1];
        }
        else
        {
            w[1] = ab[*kd + 1 + ab_dim1];
        }
        if(wantz)
        {
            z__[z_dim1 + 1] = 1.f;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants. */
    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    /* Scale matrix to allowable range, if necessary. */
    anrm = slansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1]);
    iscale = 0;
    if(anrm > 0.f && anrm < rmin)
    {
        iscale = 1;
        sigma = rmin / anrm;
    }
    else if(anrm > rmax)
    {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if(iscale == 1)
    {
        if(lower)
        {
            slascl_("B", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        else
        {
            slascl_("Q", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
    }
    /* Call SSBTRD to reduce symmetric band matrix to tridiagonal form. */
    inde = 1;
    indwrk = inde + *n;
    ssbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[inde], &z__[z_offset], ldz,
            &work[indwrk], &iinfo);
    /* For eigenvalues only, call SSTERF. For eigenvectors, call SSTEQR. */
    if(!wantz)
    {
        ssterf_(n, &w[1], &work[inde], info);
    }
    else
    {
        ssteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[indwrk], info);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if(iscale == 1)
    {
        if(*info == 0)
        {
            imax = *n;
        }
        else
        {
            imax = *info - 1;
        }
        r__1 = 1.f / sigma;
        sscal_(&imax, &r__1, &w[1], &c__1);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSBEV */
}
/* ssbev_ */
