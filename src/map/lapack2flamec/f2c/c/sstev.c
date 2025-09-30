/* ../netlib/sstev.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief <b> SSTEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for OTHER m atrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSTEV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstev.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstev.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstev.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ */
/* INTEGER INFO, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEV computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric tridiagonal matrix A. */
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
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the n diagonal elements of the tridiagonal matrix */
/* > A. */
/* > On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* > matrix A, stored in elements 1 to N-1 of E. */
/* > On exit, the contents of E are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* > eigenvectors of the matrix A, with the i-th column of Z */
/* > holding the eigenvector associated with D(i). */
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
/* > WORK is REAL array, dimension (fla_max(1,2*N-2)) */
/* > If JOBZ = 'N', WORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the algorithm failed to converge;
i */
/* > off-diagonal elements of E did not converge to zero. */
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
/** Generated wrapper function */
void sstev_(char *jobz, aocl_int_t *n, real *d__, real *e, real *z__, aocl_int_t *ldz, real *work,
            aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sstev(jobz, n, d__, e, z__, ldz, work, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sstev(jobz, &n_64, d__, e, z__, &ldz_64, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sstev(char *jobz, aocl_int64_t *n, real *d__, real *e, real *z__,
                       aocl_int64_t *ldz, real *work, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sstev inputs: jobz %c, n %" FLA_IS ", ldz %" FLA_IS "", *jobz, *n, *ldz);
    /* System generated locals */
    aocl_int64_t z_dim1, z_offset, i__1;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real eps;
    aocl_int64_t imax;
    real rmin, rmax, tnrm, sigma;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical wantz;
    aocl_int64_t iscale;
    extern real slamch_(char *);
    real safmin;
    real bignum;
    real smlnum;
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
    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SSTEV ", &i__1, (ftnlen)6);
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
    iscale = 0;
    tnrm = aocl_lapack_slanst("M", n, &d__[1], &e[1]);
    if(tnrm > 0.f && tnrm < rmin)
    {
        iscale = 1;
        sigma = rmin / tnrm;
    }
    else if(tnrm > rmax)
    {
        iscale = 1;
        sigma = rmax / tnrm;
    }
    if(iscale == 1)
    {
        aocl_blas_sscal(n, &sigma, &d__[1], &c__1);
        i__1 = *n - 1;
        aocl_blas_sscal(&i__1, &sigma, &e[1], &c__1);
    }
    /* For eigenvalues only, call SSTERF. For eigenvalues and */
    /* eigenvectors, call SSTEQR. */
    if(!wantz)
    {
        aocl_lapack_ssterf(n, &d__[1], &e[1], info);
    }
    else
    {
        aocl_lapack_ssteqr("I", n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], info);
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
        aocl_blas_sscal(&imax, &r__1, &d__[1], &c__1);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSTEV */
}
/* sstev_ */
