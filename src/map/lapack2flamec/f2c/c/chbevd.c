/* ./chbevd.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {0.f, 0.f};
static scomplex c_b2 = {1.f, 0.f};
static real c_b13 = 1.f;
static aocl_int64_t c__1 = 1;
/* > \brief <b> CHBEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHBEVD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/* LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEVD computes all the eigenvalues and, optionally, eigenvectors of */
/* > a scomplex Hermitian band matrix A. If eigenvectors are desired, it */
/* > uses a divide and conquer algorithm. */
/* > */
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
/* > AB is COMPLEX array, dimension (LDAB, N) */
/* > On entry, the upper or lower triangle of the Hermitian band */
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
/* > Z is COMPLEX array, dimension (LDZ, N) */
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
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If N <= 1, LWORK must be at least 1. */
/* > If JOBZ = 'N' and N > 1, LWORK must be at least N. */
/* > If JOBZ = 'V' and N > 1, LWORK must be at least 2*N**2. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal sizes of the WORK, RWORK and */
/* > IWORK arrays, returns these values as the first entries of */
/* > the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, */
/* > dimension (LRWORK) */
/* > On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of array RWORK. */
/* > If N <= 1, LRWORK must be at least 1. */
/* > If JOBZ = 'N' and N > 1, LRWORK must be at least N. */
/* > If JOBZ = 'V' and N > 1, LRWORK must be at least */
/* > 1 + 5*N + 2*N**2. */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of array IWORK. */
/* > If JOBZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* > If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N . */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
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
/* > \ingroup hbevd */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chbevd_(char *jobz, char *uplo, aocl_int_t *n, aocl_int_t *kd, scomplex *ab, aocl_int_t *ldab,
             real *w, scomplex *z__, aocl_int_t *ldz, scomplex *work, aocl_int_t *lwork, real *rwork,
             aocl_int_t *lrwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chbevd(jobz, uplo, n, kd, ab, ldab, w, z__, ldz, work, lwork, rwork, lrwork, iwork,
                       liwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t kd_64 = *kd;
    aocl_int64_t ldab_64 = *ldab;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t lrwork_64 = *lrwork;
    aocl_int64_t liwork_64 = *liwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chbevd(jobz, uplo, &n_64, &kd_64, ab, &ldab_64, w, z__, &ldz_64, work, &lwork_64,
                       rwork, &lrwork_64, iwork, &liwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_chbevd(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab,
                        aocl_int64_t *ldab, real *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work,
                        aocl_int64_t *lwork, real *rwork, aocl_int64_t *lrwork, aocl_int_t *iwork,
                        aocl_int64_t *liwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "chbevd inputs: jobz %c, uplo %c, n %lld, kd %lld, ldab %lld, ldz %lld, lwork %lld, "
             "lrwork %lld, liwork %lld",
             *jobz, *uplo, *n, *kd, *ldab, *ldz, *lwork, *lrwork, *liwork);
#else
    snprintf(buffer, 256,
             "chbevd inputs: jobz %c, uplo %c, n %d, kd %d, ldab %d, ldz %d, lwork %d, lrwork %d, "
             "liwork %d",
             *jobz, *uplo, *n, *kd, *ldab, *ldz, *lwork, *lrwork, *liwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real eps;
    aocl_int64_t inde;
    real anrm;
    aocl_int64_t imax;
    real rmin, rmax;
    aocl_int64_t llwk2;
    real sigma;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo;
    aocl_int64_t lwmin;
    logical lower;
    aocl_int64_t llrwk;
    logical wantz;
    aocl_int64_t indwk2;
    aocl_int64_t iscale;
    extern real slamch_(char *);
    real safmin;
    real bignum;
    aocl_int64_t indwrk, liwmin;
    aocl_int64_t lrwmin;
    real smlnum;
    logical lquery;
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    --rwork;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    lquery = *lwork == -1 || *liwork == -1 || *lrwork == -1;
    *info = 0;
    if(*n <= 1)
    {
        lwmin = 1;
        lrwmin = 1;
        liwmin = 1;
    }
    else
    {
        if(wantz)
        {
            /* Computing 2nd power */
            i__1 = *n;
            lwmin = i__1 * i__1 << 1;
            /* Computing 2nd power */
            i__1 = *n;
            lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
            liwmin = *n * 5 + 3;
        }
        else
        {
            lwmin = *n;
            lrwmin = *n;
            liwmin = 1;
        }
    }
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
    if(*info == 0)
    {
        r__1 = aocl_lapack_sroundup_lwork(&lwmin);
        work[1].real = r__1;
        work[1].imag = 0.f; // , expr subst
        rwork[1] = (real)lrwmin;
        iwork[1] = (aocl_int_t)(liwmin);
        if(*lwork < lwmin && !lquery)
        {
            *info = -11;
        }
        else if(*lrwork < lrwmin && !lquery)
        {
            *info = -13;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -15;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CHBEVD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*n == 1)
    {
        i__1 = ab_dim1 + 1;
        w[1] = ab[i__1].real;
        if(wantz)
        {
            i__1 = z_dim1 + 1;
            z__[i__1].real = 1.f;
            z__[i__1].imag = 0.f; // , expr subst
        }
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
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
    anrm = aocl_lapack_clanhb("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1]);
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
            aocl_lapack_clascl("B", kd, kd, &c_b13, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        else
        {
            aocl_lapack_clascl("Q", kd, kd, &c_b13, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
    }
    /* Call CHBTRD to reduce Hermitian band matrix to tridiagonal form. */
    inde = 1;
    indwrk = inde + *n;
    indwk2 = *n * *n + 1;
    llwk2 = *lwork - indwk2 + 1;
    llrwk = *lrwork - indwrk + 1;
    aocl_lapack_chbtrd(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &rwork[inde], &z__[z_offset],
                       ldz, &work[1], &iinfo);
    /* For eigenvalues only, call SSTERF. For eigenvectors, call CSTEDC. */
    if(!wantz)
    {
        aocl_lapack_ssterf(n, &w[1], &rwork[inde], info);
    }
    else
    {
        aocl_lapack_cstedc("I", n, &w[1], &rwork[inde], &work[1], n, &work[indwk2], &llwk2,
                           &rwork[indwrk], &llrwk, &iwork[1], liwork, info);
        aocl_blas_cgemm("N", "N", n, n, n, &c_b2, &z__[z_offset], ldz, &work[1], n, &c_b1,
                        &work[indwk2], n);
        aocl_lapack_clacpy("A", n, n, &work[indwk2], n, &z__[z_offset], ldz);
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
        aocl_blas_sscal(&imax, &r__1, &w[1], &c__1);
    }
    r__1 = aocl_lapack_sroundup_lwork(&lwmin);
    work[1].real = r__1;
    work[1].imag = 0.f; // , expr subst
    rwork[1] = (real)lrwmin;
    iwork[1] = (aocl_int_t)(liwmin);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHBEVD */
}
/* chbevd_ */
