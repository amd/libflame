/* ./strsna.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
/* > \brief \b STRSNA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STRSNA + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsna.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsna.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsna.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/* LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER HOWMNY, JOB */
/* INTEGER INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL SELECT( * ) */
/* INTEGER IWORK( * ) */
/* REAL S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ), */
/* $ VR( LDVR, * ), WORK( LDWORK, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRSNA estimates reciprocal condition numbers for specified */
/* > eigenvalues and/or right eigenvectors of a real upper */
/* > quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q */
/* > orthogonal). */
/* > */
/* > T must be in Schur canonical form (as returned by SHSEQR), that is, */
/* > block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
each */
/* > 2-by-2 diagonal block has its diagonal elements equal and its */
/* > off-diagonal elements of opposite sign. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is CHARACTER*1 */
/* > Specifies whether condition numbers are required for */
/* > eigenvalues (S) or eigenvectors (SEP): */
/* > = 'E': for eigenvalues only (S);
 */
/* > = 'V': for eigenvectors only (SEP);
 */
/* > = 'B': for both eigenvalues and eigenvectors (S and SEP). */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* > HOWMNY is CHARACTER*1 */
/* > = 'A': compute condition numbers for all eigenpairs;
 */
/* > = 'S': compute condition numbers for selected eigenpairs */
/* > specified by the array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* > SELECT is LOGICAL array, dimension (N) */
/* > If HOWMNY = 'S', SELECT specifies the eigenpairs for which */
/* > condition numbers are required. To select condition numbers */
/* > for the eigenpair corresponding to a real eigenvalue w(j), */
/* > SELECT(j) must be set to .TRUE.. To select condition numbers */
/* > corresponding to a complex conjugate pair of eigenvalues w(j) */
/* > and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be */
/* > set to .TRUE.. */
/* > If HOWMNY = 'A', SELECT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is REAL array, dimension (LDT,N) */
/* > The upper quasi-triangular matrix T, in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* > VL is REAL array, dimension (LDVL,M) */
/* > If JOB = 'E' or 'B', VL must contain left eigenvectors of T */
/* > (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
/* > eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* > must be stored in consecutive columns of VL, as returned by */
/* > SHSEIN or STREVC. */
/* > If JOB = 'V', VL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the array VL. */
/* > LDVL >= 1;
and if JOB = 'E' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] VR */
/* > \verbatim */
/* > VR is REAL array, dimension (LDVR,M) */
/* > If JOB = 'E' or 'B', VR must contain right eigenvectors of T */
/* > (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
/* > eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* > must be stored in consecutive columns of VR, as returned by */
/* > SHSEIN or STREVC. */
/* > If JOB = 'V', VR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the array VR. */
/* > LDVR >= 1;
and if JOB = 'E' or 'B', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is REAL array, dimension (MM) */
/* > If JOB = 'E' or 'B', the reciprocal condition numbers of the */
/* > selected eigenvalues, stored in consecutive elements of the */
/* > array. For a complex conjugate pair of eigenvalues two */
/* > consecutive elements of S are set to the same value. Thus */
/* > S(j), SEP(j), and the j-th columns of VL and VR all */
/* > correspond to the same eigenpair (but not in general the */
/* > j-th eigenpair, unless all eigenpairs are selected). */
/* > If JOB = 'V', S is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SEP */
/* > \verbatim */
/* > SEP is REAL array, dimension (MM) */
/* > If JOB = 'V' or 'B', the estimated reciprocal condition */
/* > numbers of the selected eigenvectors, stored in consecutive */
/* > elements of the array. For a complex eigenvector two */
/* > consecutive elements of SEP are set to the same value. If */
/* > the eigenvalues cannot be reordered to compute SEP(j), SEP(j) */
/* > is set to 0;
this can only occur when the true value would be */
/* > very small anyway. */
/* > If JOB = 'E', SEP is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* > MM is INTEGER */
/* > The number of elements in the arrays S (if JOB = 'E' or 'B') */
/* > and/or SEP (if JOB = 'V' or 'B'). MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of elements of the arrays S and/or SEP actually */
/* > used to store the estimated condition numbers. */
/* > If HOWMNY = 'A', M is set to N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (LDWORK,N+6) */
/* > If JOB = 'E', WORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* > LDWORK is INTEGER */
/* > The leading dimension of the array WORK. */
/* > LDWORK >= 1;
and if JOB = 'V' or 'B', LDWORK >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (2*(N-1)) */
/* > If JOB = 'E', IWORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup trsna */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The reciprocal of the condition number of an eigenvalue lambda is */
/* > defined as */
/* > */
/* > S(lambda) = |v**T*u| / (norm(u)*norm(v)) */
/* > */
/* > where u and v are the right and left eigenvectors of T corresponding */
/* > to lambda;
v**T denotes the transpose of v, and norm(u) */
/* > denotes the Euclidean norm. These reciprocal condition numbers always */
/* > lie between zero (very badly conditioned) and one (very well */
/* > conditioned). If n = 1, S(lambda) is defined to be 1. */
/* > */
/* > An approximate error bound for a computed eigenvalue W(i) is given by */
/* > */
/* > EPS * norm(T) / S(i) */
/* > */
/* > where EPS is the machine precision. */
/* > */
/* > The reciprocal of the condition number of the right eigenvector u */
/* > corresponding to lambda is defined as follows. Suppose */
/* > */
/* > T = ( lambda c ) */
/* > ( 0 T22 ) */
/* > */
/* > Then the reciprocal condition number is */
/* > */
/* > SEP( lambda, T22 ) = sigma-min( T22 - lambda*I ) */
/* > */
/* > where sigma-min denotes the smallest singular value. We approximate */
/* > the smallest singular value by the reciprocal of an estimate of the */
/* > one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is */
/* > defined to be f2c_abs(T(1,1)). */
/* > */
/* > An approximate error bound for a computed right eigenvector VR(i) */
/* > is given by */
/* > */
/* > EPS * norm(T) / SEP(i) */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void strsna_(char *job, char *howmny, logical *select, integer *n, real *t, integer *ldt, real *vl,
             integer *ldvl, real *vr, integer *ldvr, real *s, real *sep, integer *mm, integer *m,
             real *work, integer *ldwork, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
             "strsna inputs: job %c, howmny %c, n %" FLA_IS ", ldt %" FLA_IS ", ldvl %" FLA_IS
             ", ldvr %" FLA_IS ", mm %" FLA_IS ", ldwork %" FLA_IS "",
             *job, *howmny, *n, *ldt, *ldvl, *ldvr, *mm, *ldwork);
    /* System generated locals */
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, work_dim1, work_offset, i__1,
        i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, n2;
    real cs;
    integer nn, ks;
    real sn, mu, eps, est;
    integer kase;
    real cond;
    logical pair;
    integer ierr;
    real dumm, prod;
    integer ifst;
    real lnrm;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    integer ilst;
    real rnrm, prod1, prod2;
    extern real snrm2_(integer *, real *, integer *);
    real scale, delta;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    logical wants;
    real dummy[1];
    extern /* Subroutine */
        void
        slacn2_(integer *, real *, real *, integer *, real *, integer *, integer *);
    extern real slapy2_(real *, real *), slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    logical wantbh;
    extern /* Subroutine */
        void
        slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    logical somcon;
    extern /* Subroutine */
        void
        slaqtr_(logical *, logical *, integer *, real *, integer *, real *, real *, real *, real *,
                real *, integer *),
        strexc_(char *, integer *, real *, integer *, real *, integer *, integer *, integer *,
                real *, integer *);
    real smlnum;
    logical wantsp;
    /* -- LAPACK computational routine -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode and test the input parameters */
    /* Parameter adjustments */
    --select;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --s;
    --sep;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    --iwork;
    /* Function Body */
    wantbh = lsame_(job, "B", 1, 1);
    wants = lsame_(job, "E", 1, 1) || wantbh;
    wantsp = lsame_(job, "V", 1, 1) || wantbh;
    somcon = lsame_(howmny, "S", 1, 1);
    *info = 0;
    dumm = 0.f;
    if(!wants && !wantsp)
    {
        *info = -1;
    }
    else if(!lsame_(howmny, "A", 1, 1) && !somcon)
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*ldt < fla_max(1, *n))
    {
        *info = -6;
    }
    else if(*ldvl < 1 || wants && *ldvl < *n)
    {
        *info = -8;
    }
    else if(*ldvr < 1 || wants && *ldvr < *n)
    {
        *info = -10;
    }
    else
    {
        /* Set M to the number of eigenpairs for which condition numbers */
        /* are required, and test MM. */
        if(somcon)
        {
            *m = 0;
            pair = FALSE_;
            i__1 = *n;
            for(k = 1; k <= i__1; ++k)
            {
                if(pair)
                {
                    pair = FALSE_;
                }
                else
                {
                    if(k < *n)
                    {
                        if(t[k + 1 + k * t_dim1] == 0.f)
                        {
                            if(select[k])
                            {
                                ++(*m);
                            }
                        }
                        else
                        {
                            pair = TRUE_;
                            if(select[k] || select[k + 1])
                            {
                                *m += 2;
                            }
                        }
                    }
                    else
                    {
                        if(select[*n])
                        {
                            ++(*m);
                        }
                    }
                }
                /* L10: */
            }
        }
        else
        {
            *m = *n;
        }
        if(*mm < *m)
        {
            *info = -13;
        }
        else if(*ldwork < 1 || wantsp && *ldwork < *n)
        {
            *info = -16;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STRSNA", &i__1, (ftnlen)6);
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
        if(somcon)
        {
            if(!select[1])
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
        }
        if(wants)
        {
            s[1] = 1.f;
        }
        if(wantsp)
        {
            sep[1] = (r__1 = t[t_dim1 + 1], f2c_abs(r__1));
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants */
    eps = slamch_("P");
    smlnum = slamch_("S") / eps;
    bignum = 1.f / smlnum;
    ks = 0;
    pair = FALSE_;
    i__1 = *n;
    for(k = 1; k <= i__1; ++k)
    {
        /* Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block. */
        if(pair)
        {
            pair = FALSE_;
            goto L60;
        }
        else
        {
            if(k < *n)
            {
                pair = t[k + 1 + k * t_dim1] != 0.f;
            }
        }
        /* Determine whether condition numbers are required for the k-th */
        /* eigenpair. */
        if(somcon)
        {
            if(pair)
            {
                if(!select[k] && !select[k + 1])
                {
                    goto L60;
                }
            }
            else
            {
                if(!select[k])
                {
                    goto L60;
                }
            }
        }
        ++ks;
        if(wants)
        {
            /* Compute the reciprocal condition number of the k-th */
            /* eigenvalue. */
            if(!pair)
            {
                /* Real eigenvalue. */
                prod = sdot_(n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
                rnrm = snrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
                lnrm = snrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
                s[ks] = f2c_abs(prod) / (rnrm * lnrm);
            }
            else
            {
                /* Complex eigenvalue. */
                prod1 = sdot_(n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
                prod1 += sdot_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
                               &c__1);
                prod2 = sdot_(n, &vl[ks * vl_dim1 + 1], &c__1, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
                prod2 -= sdot_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1, &vr[ks * vr_dim1 + 1], &c__1);
                r__1 = snrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
                r__2 = snrm2_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
                rnrm = slapy2_(&r__1, &r__2);
                r__1 = snrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
                r__2 = snrm2_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1);
                lnrm = slapy2_(&r__1, &r__2);
                cond = slapy2_(&prod1, &prod2) / (rnrm * lnrm);
                s[ks] = cond;
                s[ks + 1] = cond;
            }
        }
        if(wantsp)
        {
            /* Estimate the reciprocal condition number of the k-th */
            /* eigenvector. */
            /* Copy the matrix T to the array WORK and swap the diagonal */
            /* block beginning at T(k,k) to the (1,1) position. */
            slacpy_("Full", n, n, &t[t_offset], ldt, &work[work_offset], ldwork);
            ifst = k;
            ilst = 1;
            strexc_("No Q", n, &work[work_offset], ldwork, dummy, &c__1, &ifst, &ilst,
                    &work[(*n + 1) * work_dim1 + 1], &ierr);
            if(ierr == 1 || ierr == 2)
            {
                /* Could not swap because blocks not well separated */
                scale = 1.f;
                est = bignum;
            }
            else
            {
                /* Reordering successful */
                if(work[work_dim1 + 2] == 0.f)
                {
                    /* Form C = T22 - lambda*I in WORK(2:N,2:N). */
                    i__2 = *n;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        work[i__ + i__ * work_dim1] -= work[work_dim1 + 1];
                        /* L20: */
                    }
                    n2 = 1;
                    nn = *n - 1;
                }
                else
                {
                    /* Triangularize the 2 by 2 block by unitary */
                    /* transformation U = [ cs i*ss ] */
                    /* [ i*ss cs ]. */
                    /* such that the (1,1) position of WORK is complex */
                    /* eigenvalue lambda with positive imaginary part. (2,2) */
                    /* position of WORK is the complex eigenvalue lambda */
                    /* with negative imaginary part. */
                    mu = sqrt((r__1 = work[(work_dim1 << 1) + 1], f2c_abs(r__1)))
                         * sqrt((r__2 = work[work_dim1 + 2], f2c_abs(r__2)));
                    delta = slapy2_(&mu, &work[work_dim1 + 2]);
                    cs = mu / delta;
                    sn = -work[work_dim1 + 2] / delta;
                    /* Form */
                    /* C**T = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ] */
                    /* [ mu ] */
                    /* [ .. ] */
                    /* [ .. ] */
                    /* [ mu ] */
                    /* where C**T is transpose of matrix C, */
                    /* and RWORK is stored starting in the N+1-st column of */
                    /* WORK. */
                    i__2 = *n;
                    for(j = 3; j <= i__2; ++j)
                    {
                        work[j * work_dim1 + 2] = cs * work[j * work_dim1 + 2];
                        work[j + j * work_dim1] -= work[work_dim1 + 1];
                        /* L30: */
                    }
                    work[(work_dim1 << 1) + 2] = 0.f;
                    work[(*n + 1) * work_dim1 + 1] = mu * 2.f;
                    i__2 = *n - 1;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        work[i__ + (*n + 1) * work_dim1] = sn * work[(i__ + 1) * work_dim1 + 1];
                        /* L40: */
                    }
                    n2 = 2;
                    nn = *n - 1 << 1;
                }
                /* Estimate norm(inv(C**T)) */
                est = 0.f;
                kase = 0;
            L50:
                slacn2_(&nn, &work[(*n + 2) * work_dim1 + 1], &work[(*n + 4) * work_dim1 + 1],
                        &iwork[1], &est, &kase, isave);
                if(kase != 0)
                {
                    if(kase == 1)
                    {
                        if(n2 == 1)
                        {
                            /* Real eigenvalue: solve C**T*x = scale*c. */
                            i__2 = *n - 1;
                            slaqtr_(&c_true, &c_true, &i__2, &work[(work_dim1 << 1) + 2], ldwork,
                                    dummy, &dumm, &scale, &work[(*n + 4) * work_dim1 + 1],
                                    &work[(*n + 6) * work_dim1 + 1], &ierr);
                        }
                        else
                        {
                            /* Complex eigenvalue: solve */
                            /* C**T*(p+iq) = scale*(c+id) in real arithmetic. */
                            i__2 = *n - 1;
                            slaqtr_(&c_true, &c_false, &i__2, &work[(work_dim1 << 1) + 2], ldwork,
                                    &work[(*n + 1) * work_dim1 + 1], &mu, &scale,
                                    &work[(*n + 4) * work_dim1 + 1],
                                    &work[(*n + 6) * work_dim1 + 1], &ierr);
                        }
                    }
                    else
                    {
                        if(n2 == 1)
                        {
                            /* Real eigenvalue: solve C*x = scale*c. */
                            i__2 = *n - 1;
                            slaqtr_(&c_false, &c_true, &i__2, &work[(work_dim1 << 1) + 2], ldwork,
                                    dummy, &dumm, &scale, &work[(*n + 4) * work_dim1 + 1],
                                    &work[(*n + 6) * work_dim1 + 1], &ierr);
                        }
                        else
                        {
                            /* Complex eigenvalue: solve */
                            /* C*(p+iq) = scale*(c+id) in real arithmetic. */
                            i__2 = *n - 1;
                            slaqtr_(&c_false, &c_false, &i__2, &work[(work_dim1 << 1) + 2], ldwork,
                                    &work[(*n + 1) * work_dim1 + 1], &mu, &scale,
                                    &work[(*n + 4) * work_dim1 + 1],
                                    &work[(*n + 6) * work_dim1 + 1], &ierr);
                        }
                    }
                    goto L50;
                }
            }
            sep[ks] = scale / fla_max(est, smlnum);
            if(pair)
            {
                sep[ks + 1] = sep[ks];
            }
        }
        if(pair)
        {
            ++ks;
        }
    L60:;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of STRSNA */
}
/* strsna_ */
