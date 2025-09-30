/* ./cgees.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c__0 = 0;
static aocl_int64_t c_n1 = -1;
/* > \brief <b> CGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
 * vectors f or GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEES + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgees.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgees.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgees.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, */
/* LDVS, WORK, LWORK, RWORK, BWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVS, SORT */
/* INTEGER INFO, LDA, LDVS, LWORK, N, SDIM */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL BWORK( * ) */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * ) */
/* .. */
/* .. Function Arguments .. */
/* LOGICAL SELECT */
/* EXTERNAL SELECT */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEES computes for an N-by-N scomplex nonsymmetric matrix A, the */
/* > eigenvalues, the Schur form T, and, optionally, the matrix of Schur */
/* > vectors Z. This gives the Schur factorization A = Z*T*(Z**H). */
/* > */
/* > Optionally, it also orders the eigenvalues on the diagonal of the */
/* > Schur form so that selected eigenvalues are at the top left. */
/* > The leading columns of Z then form an orthonormal basis for the */
/* > invariant subspace corresponding to the selected eigenvalues. */
/* > */
/* > A scomplex matrix is in Schur form if it is upper triangular. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVS */
/* > \verbatim */
/* > JOBVS is CHARACTER*1 */
/* > = 'N': Schur vectors are not computed;
 */
/* > = 'V': Schur vectors are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] SORT */
/* > \verbatim */
/* > SORT is CHARACTER*1 */
/* > Specifies whether or not to order the eigenvalues on the */
/* > diagonal of the Schur form. */
/* > = 'N': Eigenvalues are not ordered: */
/* > = 'S': Eigenvalues are ordered (see SELECT). */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* > SELECT is a LOGICAL FUNCTION of one COMPLEX argument */
/* > SELECT must be declared EXTERNAL in the calling subroutine. */
/* > If SORT = 'S', SELECT is used to select eigenvalues to order */
/* > to the top left of the Schur form. */
/* > IF SORT = 'N', SELECT is not referenced. */
/* > The eigenvalue W(j) is selected if SELECT(W(j)) is true. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. */
/* > On exit, A has been overwritten by its Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SDIM */
/* > \verbatim */
/* > SDIM is INTEGER */
/* > If SORT = 'N', SDIM = 0. */
/* > If SORT = 'S', SDIM = number of eigenvalues for which */
/* > SELECT is true. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (N) */
/* > W contains the computed eigenvalues, in the same order that */
/* > they appear on the diagonal of the output Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[out] VS */
/* > \verbatim */
/* > VS is COMPLEX array, dimension (LDVS,N) */
/* > If JOBVS = 'V', VS contains the unitary matrix Z of Schur */
/* > vectors. */
/* > If JOBVS = 'N', VS is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVS */
/* > \verbatim */
/* > LDVS is INTEGER */
/* > The leading dimension of the array VS. LDVS >= 1;
if */
/* > JOBVS = 'V', LDVS >= N. */
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
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BWORK */
/* > \verbatim */
/* > BWORK is LOGICAL array, dimension (N) */
/* > Not referenced if SORT = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = i, and i is */
/* > <= N: the QR algorithm failed to compute all the */
/* > eigenvalues;
elements 1:ILO-1 and i+1:N of W */
/* > contain those eigenvalues which have converged;
 */
/* > if JOBVS = 'V', VS contains the matrix which */
/* > reduces A to its partially converged Schur form. */
/* > = N+1: the eigenvalues could not be reordered because */
/* > some eigenvalues were too close to separate (the */
/* > problem is very ill-conditioned);
 */
/* > = N+2: after reordering, roundoff changed values of */
/* > some scomplex eigenvalues so that leading */
/* > eigenvalues in the Schur form no longer satisfy */
/* > SELECT = .TRUE.. This could also be caused by */
/* > underflow due to scaling. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gees */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgees_(char *jobvs, char *sort, L_fp1 select, aocl_int_t *n, scomplex *a, aocl_int_t *lda,
            aocl_int_t *sdim, scomplex *w, scomplex *vs, aocl_int_t *ldvs, scomplex *work,
            aocl_int_t *lwork, real *rwork, logical *bwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgees(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork,
                      info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t sdim_64 = *sdim;
    aocl_int64_t ldvs_64 = *ldvs;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgees(jobvs, sort, select, &n_64, a, &lda_64, &sdim_64, w, vs, &ldvs_64, work,
                      &lwork_64, rwork, bwork, &info_64);

    *sdim = (aocl_int_t)sdim_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgees(char *jobvs, char *sort, L_fp1 select, aocl_int64_t *n, scomplex *a,
                       aocl_int64_t *lda, aocl_int64_t *sdim, scomplex *w, scomplex *vs,
                       aocl_int64_t *ldvs, scomplex *work, aocl_int64_t *lwork, real *rwork,
                       logical *bwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgees inputs: jobvs %c, sort %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldvs %" FLA_IS "",
                      *jobvs, *sort, *n, *lda, *ldvs);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    aocl_int64_t i__;
    real s;
    aocl_int64_t ihi, ilo;
    real dum[1], eps, sep;
    aocl_int64_t ibal;
    real anrm;
    aocl_int64_t ierr, itau, iwrk, icond, ieval;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical scalea;
    real cscale;
    extern real slamch_(char *);
    real bignum;
    aocl_int64_t minwrk, maxwrk;
    real smlnum;
    aocl_int64_t hswork;
    logical wantst, lquery, wantvs;
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Function Arguments .. */
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
    vs_dim1 = *ldvs;
    vs_offset = 1 + vs_dim1;
    vs -= vs_offset;
    --work;
    --rwork;
    --bwork;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    wantvs = lsame_(jobvs, "V", 1, 1);
    wantst = lsame_(sort, "S", 1, 1);
    if(!wantvs && !lsame_(jobvs, "N", 1, 1))
    {
        *info = -1;
    }
    else if(!wantst && !lsame_(sort, "N", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -6;
    }
    else if(*ldvs < 1 || wantvs && *ldvs < *n)
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
    /* HSWORK refers to the workspace preferred by CHSEQR, as */
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
            maxwrk = *n + *n * aocl_lapack_ilaenv(&c__1, "CGEHRD", " ", n, &c__1, n, &c__0);
            minwrk = *n << 1;
            aocl_lapack_chseqr("S", jobvs, n, &c__1, n, &a[a_offset], lda, &w[1], &vs[vs_offset],
                               ldvs, &work[1], &c_n1, &ieval);
            hswork = (integer)work[1].r;
            if(!wantvs)
            {
                maxwrk = fla_max(maxwrk, hswork);
            }
            else
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n
                       + (*n - 1)
                             * aocl_lapack_ilaenv(&c__1, "CUNGHR", " ", n, &c__1, n,
                                                  &c_n1); // , expr subst
                maxwrk = fla_max(i__1, i__2);
                maxwrk = fla_max(maxwrk, hswork);
            }
        }
        r__1 = aocl_lapack_sroundup_lwork(&maxwrk);
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
        if(*lwork < minwrk && !lquery)
        {
            *info = -12;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGEES ", &i__1, (ftnlen)6);
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
        *sdim = 0;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants */
    eps = slamch_("P");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = aocl_lapack_clange("M", n, n, &a[a_offset], lda, dum);
    scalea = FALSE_;
    if(anrm > 0.f && anrm < smlnum)
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
        aocl_lapack_clascl("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &ierr);
    }
    /* Permute the matrix to make it more nearly triangular */
    /* (CWorkspace: none) */
    /* (RWorkspace: need N) */
    ibal = 1;
    aocl_lapack_cgebal("P", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr);
    /* Reduce to upper Hessenberg form */
    /* (CWorkspace: need 2*N, prefer N+N*NB) */
    /* (RWorkspace: none) */
    itau = 1;
    iwrk = *n + itau;
    i__1 = *lwork - iwrk + 1;
    aocl_lapack_cgehrd(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);
    if(wantvs)
    {
        /* Copy Householder vectors to VS */
        aocl_lapack_clacpy("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs);
        /* Generate unitary matrix in VS */
        /* (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_cunghr(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk], &i__1,
                           &ierr);
    }
    *sdim = 0;
    /* Perform QR iteration, accumulating Schur vectors in VS if desired */
    /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
    /* (RWorkspace: none) */
    iwrk = itau;
    i__1 = *lwork - iwrk + 1;
    aocl_lapack_chseqr("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vs[vs_offset], ldvs,
                       &work[iwrk], &i__1, &ieval);
    if(ieval > 0)
    {
        *info = ieval;
    }
    /* Sort eigenvalues if desired */
    if(wantst && *info == 0)
    {
        if(scalea)
        {
            aocl_lapack_clascl("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &w[1], n, &ierr);
        }
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            bwork[i__] = (*select)(&w[i__]);
            /* L10: */
        }
        /* Reorder eigenvalues and transform Schur vectors */
        /* (CWorkspace: none) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_ctrsen("N", jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset], ldvs, &w[1],
                           sdim, &s, &sep, &work[iwrk], &i__1, &icond);
    }
    if(wantvs)
    {
        /* Undo balancing */
        /* (CWorkspace: none) */
        /* (RWorkspace: need N) */
        aocl_lapack_cgebak("P", "R", n, &ilo, &ihi, &rwork[ibal], n, &vs[vs_offset], ldvs, &ierr);
    }
    if(scalea)
    {
        /* Undo scaling for the Schur form of A */
        aocl_lapack_clascl("U", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, &ierr);
        i__1 = *lda + 1;
        aocl_blas_ccopy(n, &a[a_offset], &i__1, &w[1], &c__1);
    }
    r__1 = aocl_lapack_sroundup_lwork(&maxwrk);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of CGEES */
}
/* cgees_ */
