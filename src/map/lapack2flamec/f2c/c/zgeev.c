/* ./zgeev.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c__0 = 0;
static aocl_int64_t c_n1 = -1;
/* > \brief <b> ZGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for GE matr ices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGEEV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeev.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeev.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeev.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, */
/* WORK, LWORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVL, JOBVR */
/* INTEGER INFO, LDA, LDVL, LDVR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION RWORK( * ) */
/* COMPLEX*16 A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), */
/* $ W( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEEV computes for an N-by-N scomplex nonsymmetric matrix A, the */
/* > eigenvalues and, optionally, the left and/or right eigenvectors. */
/* > */
/* > The right eigenvector v(j) of A satisfies */
/* > A * v(j) = lambda(j) * v(j) */
/* > where lambda(j) is its eigenvalue. */
/* > The left eigenvector u(j) of A satisfies */
/* > u(j)**H * A = lambda(j) * u(j)**H */
/* > where u(j)**H denotes the conjugate transpose of u(j). */
/* > */
/* > The computed eigenvectors are normalized to have Euclidean norm */
/* > equal to 1 and largest component real. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVL */
/* > \verbatim */
/* > JOBVL is CHARACTER*1 */
/* > = 'N': left eigenvectors of A are not computed;
 */
/* > = 'V': left eigenvectors of are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* > JOBVR is CHARACTER*1 */
/* > = 'N': right eigenvectors of A are not computed;
 */
/* > = 'V': right eigenvectors of A are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. */
/* > On exit, A has been overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX*16 array, dimension (N) */
/* > W contains the computed eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* > VL is COMPLEX*16 array, dimension (LDVL,N) */
/* > If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* > after another in the columns of VL, in the same order */
/* > as their eigenvalues. */
/* > If JOBVL = 'N', VL is not referenced. */
/* > u(j) = VL(:,j), the j-th column of VL. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the array VL. LDVL >= 1;
if */
/* > JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* > VR is COMPLEX*16 array, dimension (LDVR,N) */
/* > If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* > after another in the columns of VR, in the same order */
/* > as their eigenvalues. */
/* > If JOBVR = 'N', VR is not referenced. */
/* > v(j) = VR(:,j), the j-th column of VR. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the array VR. LDVR >= 1;
if */
/* > JOBVR = 'V', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,2*N). */
/* > For good performance, LWORK must generally be larger. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = i, the QR algorithm failed to compute all the */
/* > eigenvalues, and no eigenvectors have been computed;
 */
/* > elements i+1:N of W contain eigenvalues which have */
/* > converged. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* @precisions fortran z -> c */
/* > \ingroup geev */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zgeev_(char *jobvl, char *jobvr, aocl_int_t *n, dcomplex *a, aocl_int_t *lda,
            dcomplex *w, dcomplex *vl, aocl_int_t *ldvl, dcomplex *vr,
            aocl_int_t *ldvr, dcomplex *work, aocl_int_t *lwork, doublereal *rwork,
            aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldvl_64 = *ldvl;
    aocl_int64_t ldvr_64 = *ldvr;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgeev(jobvl, jobvr, &n_64, a, &lda_64, w, vl, &ldvl_64, vr, &ldvr_64, work,
                      &lwork_64, rwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zgeev(char *jobvl, char *jobvr, aocl_int64_t *n, dcomplex *a,
                       aocl_int64_t *lda, dcomplex *w, dcomplex *vl, aocl_int64_t *ldvl,
                       dcomplex *vr, aocl_int64_t *ldvr, dcomplex *work,
                       aocl_int64_t *lwork, doublereal *rwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeev inputs: jobvl %c, jobvr %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldvl %" FLA_IS ", ldvr %" FLA_IS ", lwork %" FLA_IS "",
                      *jobvl, *jobvr, *n, *lda, *ldvl, *ldvr, *lwork);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    dcomplex z__1, z__2;
    /* Builtin functions */
    double sqrt(doublereal), d_imag(dcomplex *);
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__, k, ihi;
    doublereal scl;
    aocl_int64_t ilo;
    doublereal dum[1], eps;
    dcomplex tmp;
    aocl_int64_t lwork_trevc__, ibal;
    char side[1];
    doublereal anrm;
    aocl_int64_t ierr, itau, iwrk, nout;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical scalea;
    extern doublereal dlamch_(char *);
    doublereal cscale;
    logical select[1];
    doublereal bignum;
    aocl_int64_t minwrk, maxwrk;
    logical wantvl;
    doublereal smlnum;
    aocl_int64_t hswork, irwork;
    logical lquery, wantvr;
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    wantvl = lsame_(jobvl, "V", 1, 1);
    wantvr = lsame_(jobvr, "V", 1, 1);
    if(!wantvl && !lsame_(jobvl, "N", 1, 1))
    {
        *info = -1;
    }
    else if(!wantvr && !lsame_(jobvr, "N", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*ldvl < 1 || wantvl && *ldvl < *n)
    {
        *info = -8;
    }
    else if(*ldvr < 1 || wantvr && *ldvr < *n)
    {
        *info = -10;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* CWorkspace refers to scomplex workspace, and RWorkspace to real */
    /* workspace. NB refers to the optimal block size for the */
    /* immediately following subroutine, as returned by ILAENV. */
    /* HSWORK refers to the workspace preferred by ZHSEQR, as */
    /* calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
    /* the worst case.) */
    if(*info == 0)
    {
        if(*n == 0)
        {
            minwrk = 1;
            maxwrk = 1;
        }
        else
        {
            maxwrk = *n + *n * aocl_lapack_ilaenv(&c__1, "ZGEHRD", " ", n, &c__1, n, &c__0);
            minwrk = *n << 1;
            if(wantvl)
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n
                       + (*n - 1)
                             * aocl_lapack_ilaenv(&c__1, "ZUNGHR", " ", n, &c__1, n,
                                                  &c_n1); // , expr subst
                maxwrk = fla_max(i__1, i__2);
                aocl_lapack_ztrevc3("L", "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
                                    &vr[vr_offset], ldvr, n, &nout, &work[1], &c_n1, &rwork[1],
                                    &c_n1, &ierr);
                lwork_trevc__ = (integer)work[1].real;
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + lwork_trevc__; // , expr subst
                maxwrk = fla_max(i__1, i__2);
                aocl_lapack_zhseqr("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vl[vl_offset],
                                   ldvl, &work[1], &c_n1, info);
            }
            else if(wantvr)
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n
                       + (*n - 1)
                             * aocl_lapack_ilaenv(&c__1, "ZUNGHR", " ", n, &c__1, n,
                                                  &c_n1); // , expr subst
                maxwrk = fla_max(i__1, i__2);
                aocl_lapack_ztrevc3("R", "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
                                    &vr[vr_offset], ldvr, n, &nout, &work[1], &c_n1, &rwork[1],
                                    &c_n1, &ierr);
                lwork_trevc__ = (integer)work[1].real;
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + lwork_trevc__; // , expr subst
                maxwrk = fla_max(i__1, i__2);
                aocl_lapack_zhseqr("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[vr_offset],
                                   ldvr, &work[1], &c_n1, info);
            }
            else
            {
                aocl_lapack_zhseqr("E", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[vr_offset],
                                   ldvr, &work[1], &c_n1, info);
            }
            hswork = (integer)work[1].real;
            /* Computing MAX */
            i__1 = fla_max(maxwrk, hswork);
            maxwrk = fla_max(i__1, minwrk);
        }
        work[1].real = (doublereal)maxwrk;
        work[1].imag = 0.; // , expr subst
        if(*lwork < minwrk && !lquery)
        {
            *info = -12;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZGEEV ", &i__1, (ftnlen)6);
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
    /* Get machine constants */
    eps = dlamch_("P");
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = aocl_lapack_zlange("M", n, n, &a[a_offset], lda, dum);
    scalea = FALSE_;
    if(anrm > 0. && anrm < smlnum)
    {
        scalea = TRUE_;
        cscale = smlnum;
    }
    else if(anrm > bignum)
    {
        scalea = TRUE_;
        cscale = bignum;
    }
    if(scalea)
    {
        aocl_lapack_zlascl("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &ierr);
    }
    /* Balance the matrix */
    /* (CWorkspace: none) */
    /* (RWorkspace: need N) */
    ibal = 1;
    aocl_lapack_zgebal("B", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr);
    /* Reduce to upper Hessenberg form */
    /* (CWorkspace: need 2*N, prefer N+N*NB) */
    /* (RWorkspace: none) */
    itau = 1;
    iwrk = itau + *n;
    i__1 = *lwork - iwrk + 1;
    aocl_lapack_zgehrd(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);
    if(wantvl)
    {
        /* Want left eigenvectors */
        /* Copy Householder vectors to VL */
        *(unsigned char *)side = 'L';
        aocl_lapack_zlacpy("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl);
        /* Generate unitary matrix in VL */
        /* (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_zunghr(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &i__1,
                           &ierr);
        /* Perform QR iteration, accumulating Schur vectors in VL */
        /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
        /* (RWorkspace: none) */
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_zhseqr("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vl[vl_offset], ldvl,
                           &work[iwrk], &i__1, info);
        if(wantvr)
        {
            /* Want left and right eigenvectors */
            /* Copy Schur vectors to VR */
            *(unsigned char *)side = 'B';
            aocl_lapack_zlacpy("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
        }
    }
    else if(wantvr)
    {
        /* Want right eigenvectors */
        /* Copy Householder vectors to VR */
        *(unsigned char *)side = 'R';
        aocl_lapack_zlacpy("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr);
        /* Generate unitary matrix in VR */
        /* (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_zunghr(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &i__1,
                           &ierr);
        /* Perform QR iteration, accumulating Schur vectors in VR */
        /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
        /* (RWorkspace: none) */
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_zhseqr("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[vr_offset], ldvr,
                           &work[iwrk], &i__1, info);
    }
    else
    {
        /* Compute eigenvalues only */
        /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
        /* (RWorkspace: none) */
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_zhseqr("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[vr_offset], ldvr,
                           &work[iwrk], &i__1, info);
    }
    /* If INFO .NE. 0 from ZHSEQR, then quit */
    if(*info != 0)
    {
        goto L50;
    }
    if(wantvl || wantvr)
    {
        /* Compute left and/or right eigenvectors */
        /* (CWorkspace: need 2*N, prefer N + 2*N*NB) */
        /* (RWorkspace: need 2*N) */
        irwork = ibal + *n;
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_ztrevc3(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
                            &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &i__1, &rwork[irwork], n,
                            &ierr);
    }
    if(wantvl)
    {
        /* Undo balancing of left eigenvectors */
        /* (CWorkspace: none) */
        /* (RWorkspace: need N) */
        aocl_lapack_zgebak("B", "L", n, &ilo, &ihi, &rwork[ibal], n, &vl[vl_offset], ldvl, &ierr);
        /* Normalize left eigenvectors and make largest component real */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            scl = 1. / aocl_blas_dznrm2(n, &vl[i__ * vl_dim1 + 1], &c__1);
            aocl_blas_zdscal(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
            i__2 = *n;
            for(k = 1; k <= i__2; ++k)
            {
                i__3 = k + i__ * vl_dim1;
                /* Computing 2nd power */
                d__1 = vl[i__3].real;
                /* Computing 2nd power */
                d__2 = d_imag(&vl[k + i__ * vl_dim1]);
                rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
                /* L10: */
            }
            k = aocl_blas_idamax(n, &rwork[irwork], &c__1);
            d_cnjg(&z__2, &vl[k + i__ * vl_dim1]);
            d__1 = sqrt(rwork[irwork + k - 1]);
            z__1.real = z__2.real / d__1;
            z__1.imag = z__2.imag / d__1; // , expr subst
            tmp.real = z__1.real;
            tmp.imag = z__1.imag; // , expr subst
            aocl_blas_zscal(n, &tmp, &vl[i__ * vl_dim1 + 1], &c__1);
            i__2 = k + i__ * vl_dim1;
            i__3 = k + i__ * vl_dim1;
            d__1 = vl[i__3].real;
            z__1.real = d__1;
            z__1.imag = 0.; // , expr subst
            vl[i__2].real = z__1.real;
            vl[i__2].imag = z__1.imag; // , expr subst
            /* L20: */
        }
    }
    if(wantvr)
    {
        /* Undo balancing of right eigenvectors */
        /* (CWorkspace: none) */
        /* (RWorkspace: need N) */
        aocl_lapack_zgebak("B", "R", n, &ilo, &ihi, &rwork[ibal], n, &vr[vr_offset], ldvr, &ierr);
        /* Normalize right eigenvectors and make largest component real */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            scl = 1. / aocl_blas_dznrm2(n, &vr[i__ * vr_dim1 + 1], &c__1);
            aocl_blas_zdscal(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
            i__2 = *n;
            for(k = 1; k <= i__2; ++k)
            {
                i__3 = k + i__ * vr_dim1;
                /* Computing 2nd power */
                d__1 = vr[i__3].real;
                /* Computing 2nd power */
                d__2 = d_imag(&vr[k + i__ * vr_dim1]);
                rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
                /* L30: */
            }
            k = aocl_blas_idamax(n, &rwork[irwork], &c__1);
            d_cnjg(&z__2, &vr[k + i__ * vr_dim1]);
            d__1 = sqrt(rwork[irwork + k - 1]);
            z__1.real = z__2.real / d__1;
            z__1.imag = z__2.imag / d__1; // , expr subst
            tmp.real = z__1.real;
            tmp.imag = z__1.imag; // , expr subst
            aocl_blas_zscal(n, &tmp, &vr[i__ * vr_dim1 + 1], &c__1);
            i__2 = k + i__ * vr_dim1;
            i__3 = k + i__ * vr_dim1;
            d__1 = vr[i__3].real;
            z__1.real = d__1;
            z__1.imag = 0.; // , expr subst
            vr[i__2].real = z__1.real;
            vr[i__2].imag = z__1.imag; // , expr subst
            /* L40: */
        }
    }
/* Undo scaling if necessary */
L50:
    if(scalea)
    {
        i__1 = *n - *info;
        /* Computing MAX */
        i__3 = *n - *info;
        i__2 = fla_max(i__3, 1);
        aocl_lapack_zlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[*info + 1], &i__2,
                           &ierr);
        if(*info > 0)
        {
            i__1 = ilo - 1;
            aocl_lapack_zlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[1], n, &ierr);
        }
    }
    work[1].real = (doublereal)maxwrk;
    work[1].imag = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGEEV */
}
/* zgeev_ */
