#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c__0 = 0;
static doublereal c_b17 = 1.;
/* > \brief <b> dsteqr_helper computes the eigenvalues and  right eigenvectors for SY matrices</b>
 */
/* =========== DOCUMENTATION =========== */
/* ===================================================================== */
/* Subroutine */
void dsteqr_helper_(char *jobz, char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                    doublereal *w, doublereal *work, aocl_int64_t *lwork, aocl_int_t *iwork,
                    aocl_int64_t *liwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dsteqr_helper inputs: jobz %c, uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", lwork %" FLA_IS ", liwork %" FLA_IS "",
                      *jobz, *uplo, *n, *lda, *lwork, *liwork);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal eps;
    aocl_int64_t inde;
    doublereal anrm, rmin, rmax;
    aocl_int64_t lopt;
    doublereal sigma;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo, lwmin, liopt;
    logical wantz;
    aocl_int64_t indwk2, llwrk2;
    extern doublereal dlamch_(char *);
    aocl_int64_t iscale;
    doublereal safmin;
    doublereal bignum;
    aocl_int64_t indtau;
    aocl_int64_t indwrk, liwmin;
    aocl_int64_t llwork;
    doublereal smlnum;
    logical lquery;
    /* -- LAPACK driver routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --work;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    *info = 0;

    if(*info == 0)
    {
        if(*n <= 1)
        {
            liwmin = 1;
            lwmin = 1;
            lopt = lwmin;
            liopt = liwmin;
        }
        else
        {
            if(wantz)
            {
                liwmin = *n * 5 + 3;
                /* Computing 2nd power */
                i__1 = *n;
                lwmin = *n * 6 + 1 + (i__1 * i__1 << 1);
            }
            else
            {
                liwmin = 1;
                lwmin = (*n << 1) + 1;
            }
            /* Computing MAX */
            i__1 = lwmin;
            i__2 = (*n << 1)
                   + aocl_lapack_ilaenv(&c__1, "DSYEVD", uplo, n, &c_n1, &c_n1,
                                        &c_n1); // , expr subst
            lopt = fla_max(i__1, i__2);
            liopt = liwmin;
        }
        work[1] = (doublereal)lopt;
        iwork[1] = (aocl_int_t)liopt;
        if(*lwork < lwmin && !lquery)
        {
            *info = -8;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -10;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DSYEVD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
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
        w[1] = a[a_dim1 + 1];
        if(wantz)
        {
            a[a_dim1 + 1] = 1.;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants. */
    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    /* Scale matrix to allowable range, if necessary. */
    anrm = aocl_lapack_dlansy("M", uplo, n, &a[a_offset], lda, &work[1]);
    iscale = 0;
    if(anrm > 0. && anrm < rmin)
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
        aocl_lapack_dlascl(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, info);
    }
    // printf("reaching after lascl\n");
    /* Call DSYTRD to reduce symmetric matrix to tridiagonal form. */
    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;
    /* For eigenvalues only, call DSTERF. For eigenvectors, first call */
    /* DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
    /* tridiagonal matrix, then call DORMTR to multiply it by the */
    /* Householder transformations stored in A. */
    aocl_lapack_dsytrd(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &work[indwrk],
                       &llwork, &iinfo);
    lopt = (integer)((*n << 1) + work[indwrk]);
    if(!wantz)
    {
        aocl_lapack_dsterf(n, &w[1], &work[inde], info);
    }
    else
    {
        aocl_lapack_dstedc("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &llwrk2,
                           &iwork[1], liwork, info);
        aocl_lapack_dormtr("L", uplo, "N", n, n, &a[a_offset], lda, &work[indtau], &work[indwrk], n,
                           &work[indwk2], &llwrk2, &iinfo);
        aocl_lapack_dlacpy("A", n, n, &work[indwrk], n, &a[a_offset], lda);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if(iscale == 1)
    {
        d__1 = 1. / sigma;
        aocl_blas_dscal(n, &d__1, &w[1], &c__1);
    }
    work[1] = (doublereal)lopt;
    iwork[1] = (aocl_int_t)liopt;

    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DSTEQR_HELPER */
}
/* dsteqr_helper_ */
