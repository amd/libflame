/* ../netlib/dhsein.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static logical c_false = FALSE_;
static logical c_true = TRUE_;
/* > \brief \b DHSEIN */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DHSEIN + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhsein.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhsein.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhsein.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI, */
/* VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL, */
/* IFAILR, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER EIGSRC, INITV, SIDE */
/* INTEGER INFO, LDH, LDVL, LDVR, M, MM, N */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL SELECT( * ) */
/* INTEGER IFAILL( * ), IFAILR( * ) */
/* DOUBLE PRECISION H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), */
/* $ WI( * ), WORK( * ), WR( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DHSEIN uses inverse iteration to find specified right and/or left */
/* > eigenvectors of a real upper Hessenberg matrix H. */
/* > */
/* > The right eigenvector x and the left eigenvector y of the matrix H */
/* > corresponding to an eigenvalue w are defined by: */
/* > */
/* > H * x = w * x, y**h * H = w * y**h */
/* > */
/* > where y**h denotes the conjugate transpose of the vector y. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'R': compute right eigenvectors only;
 */
/* > = 'L': compute left eigenvectors only;
 */
/* > = 'B': compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] EIGSRC */
/* > \verbatim */
/* > EIGSRC is CHARACTER*1 */
/* > Specifies the source of eigenvalues supplied in (WR,WI): */
/* > = 'Q': the eigenvalues were found using DHSEQR;
thus, if */
/* > H has zero subdiagonal elements, and so is */
/* > block-triangular, then the j-th eigenvalue can be */
/* > assumed to be an eigenvalue of the block containing */
/* > the j-th row/column. This property allows DHSEIN to */
/* > perform inverse iteration on just one diagonal block. */
/* > = 'N': no assumptions are made on the correspondence */
/* > between eigenvalues and diagonal blocks. In this */
/* > case, DHSEIN must always perform inverse iteration */
/* > using the whole matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] INITV */
/* > \verbatim */
/* > INITV is CHARACTER*1 */
/* > = 'N': no initial vectors are supplied;
 */
/* > = 'U': user-supplied initial vectors are stored in the arrays */
/* > VL and/or VR. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SELECT */
/* > \verbatim */
/* > SELECT is LOGICAL array, dimension (N) */
/* > Specifies the eigenvectors to be computed. To select the */
/* > real eigenvector corresponding to a real eigenvalue WR(j), */
/* > SELECT(j) must be set to .TRUE.. To select the scomplex */
/* > eigenvector corresponding to a scomplex eigenvalue */
/* > (WR(j),WI(j)), with scomplex conjugate (WR(j+1),WI(j+1)), */
/* > either SELECT(j) or SELECT(j+1) or both must be set to */
/* > .TRUE.;
then on exit SELECT(j) is .TRUE. and SELECT(j+1) is */
/* > .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix H. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* > H is DOUBLE PRECISION array, dimension (LDH,N) */
/* > The upper Hessenberg matrix H. */
/* > If a NaN is detected in H, the routine will return with INFO=-6. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > The leading dimension of the array H. LDH >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] WR */
/* > \verbatim */
/* > WR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in] WI */
/* > \verbatim */
/* > WI is DOUBLE PRECISION array, dimension (N) */
/* > */
/* > On entry, the real and imaginary parts of the eigenvalues of */
/* > H;
a scomplex conjugate pair of eigenvalues must be stored in */
/* > consecutive elements of WR and WI. */
/* > On exit, WR may have been altered since close eigenvalues */
/* > are perturbed slightly in searching for independent */
/* > eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* > VL is DOUBLE PRECISION array, dimension (LDVL,MM) */
/* > On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must */
/* > contain starting vectors for the inverse iteration for the */
/* > left eigenvectors;
the starting vector for each eigenvector */
/* > must be in the same column(s) in which the eigenvector will */
/* > be stored. */
/* > On exit, if SIDE = 'L' or 'B', the left eigenvectors */
/* > specified by SELECT will be stored consecutively in the */
/* > columns of VL, in the same order as their eigenvalues. A */
/* > scomplex eigenvector corresponding to a scomplex eigenvalue is */
/* > stored in two consecutive columns, the first holding the real */
/* > part and the second the imaginary part. */
/* > If SIDE = 'R', VL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the array VL. */
/* > LDVL >= fla_max(1,N) if SIDE = 'L' or 'B';
LDVL >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* > VR is DOUBLE PRECISION array, dimension (LDVR,MM) */
/* > On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must */
/* > contain starting vectors for the inverse iteration for the */
/* > right eigenvectors;
the starting vector for each eigenvector */
/* > must be in the same column(s) in which the eigenvector will */
/* > be stored. */
/* > On exit, if SIDE = 'R' or 'B', the right eigenvectors */
/* > specified by SELECT will be stored consecutively in the */
/* > columns of VR, in the same order as their eigenvalues. A */
/* > scomplex eigenvector corresponding to a scomplex eigenvalue is */
/* > stored in two consecutive columns, the first holding the real */
/* > part and the second the imaginary part. */
/* > If SIDE = 'L', VR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the array VR. */
/* > LDVR >= fla_max(1,N) if SIDE = 'R' or 'B';
LDVR >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* > MM is INTEGER */
/* > The number of columns in the arrays VL and/or VR. MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of columns in the arrays VL and/or VR required to */
/* > store the eigenvectors;
each selected real eigenvector */
/* > occupies one column and each selected scomplex eigenvector */
/* > occupies two columns. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension ((N+2)*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILL */
/* > \verbatim */
/* > IFAILL is INTEGER array, dimension (MM) */
/* > If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left */
/* > eigenvector in the i-th column of VL (corresponding to the */
/* > eigenvalue w(j)) failed to converge;
IFAILL(i) = 0 if the */
/* > eigenvector converged satisfactorily. If the i-th and (i+1)th */
/* > columns of VL hold a scomplex eigenvector, then IFAILL(i) and */
/* > IFAILL(i+1) are set to the same value. */
/* > If SIDE = 'R', IFAILL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILR */
/* > \verbatim */
/* > IFAILR is INTEGER array, dimension (MM) */
/* > If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right */
/* > eigenvector in the i-th column of VR (corresponding to the */
/* > eigenvalue w(j)) failed to converge;
IFAILR(i) = 0 if the */
/* > eigenvector converged satisfactorily. If the i-th and (i+1)th */
/* > columns of VR hold a scomplex eigenvector, then IFAILR(i) and */
/* > IFAILR(i+1) are set to the same value. */
/* > If SIDE = 'L', IFAILR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, i is the number of eigenvectors which */
/* > failed to converge;
see IFAILL and IFAILR for further */
/* > details. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Each eigenvector is normalized so that the element of largest */
/* > magnitude has magnitude 1;
here the magnitude of a scomplex number */
/* > (x,y) is taken to be |x|+|y|. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dhsein_(char *side, char *eigsrc, char *initv, logical *select, aocl_int_t *n, doublereal *h__,
             aocl_int_t *ldh, doublereal *wr, doublereal *wi, doublereal *vl, aocl_int_t *ldvl,
             doublereal *vr, aocl_int_t *ldvr, aocl_int_t *mm, aocl_int_t *m, doublereal *work,
             aocl_int_t *ifaill, aocl_int_t *ifailr, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dhsein(side, eigsrc, initv, select, n, h__, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m,
                       work, ifaill, ifailr, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldh_64 = *ldh;
    aocl_int64_t ldvl_64 = *ldvl;
    aocl_int64_t ldvr_64 = *ldvr;
    aocl_int64_t mm_64 = *mm;
    aocl_int64_t m_64 = *m;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dhsein(side, eigsrc, initv, select, &n_64, h__, &ldh_64, wr, wi, vl, &ldvl_64, vr,
                       &ldvr_64, &mm_64, &m_64, work, ifaill, ifailr, &info_64);

    *m = (aocl_int_t)m_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_dhsein(char *side, char *eigsrc, char *initv, logical *select, aocl_int64_t *n,
                        doublereal *h__, aocl_int64_t *ldh, doublereal *wr, doublereal *wi,
                        doublereal *vl, aocl_int64_t *ldvl, doublereal *vr, aocl_int64_t *ldvr,
                        aocl_int64_t *mm, aocl_int64_t *m, doublereal *work, aocl_int_t *ifaill,
                        aocl_int_t *ifailr, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dhsein inputs: side %c, eigsrc %c, initv %c, n %" FLA_IS ", ldh %" FLA_IS
                      ", ldvl %" FLA_IS ", ldvr %" FLA_IS ", mm %" FLA_IS "",
                      *side, *eigsrc, *initv, *n, *ldh, *ldvl, *ldvr, *mm);
    /* System generated locals */
    aocl_int64_t h_dim1, h_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Local variables */
    aocl_int64_t i__, k, kl, kr, kln, ksi;
    doublereal wki;
    aocl_int64_t ksr;
    doublereal ulp, wkr, eps3;
    logical pair;
    doublereal unfl;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo;
    logical leftv, bothv;
    doublereal hnorm;
    extern doublereal dlamch_(char *);
    extern logical disnan_(doublereal *);
    doublereal bignum;
    logical noinit;
    aocl_int64_t ldwork;
    logical rightv, fromqr;
    doublereal smlnum;
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
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
    /* Decode and test the input parameters. */
    /* Parameter adjustments */
    --select;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --ifaill;
    --ifailr;
    /* Function Body */
    iinfo = 0;
    bothv = lsame_(side, "B", 1, 1);
    rightv = lsame_(side, "R", 1, 1) || bothv;
    leftv = lsame_(side, "L", 1, 1) || bothv;
    fromqr = lsame_(eigsrc, "Q", 1, 1);
    noinit = lsame_(initv, "N", 1, 1);
    /* Set M to the number of columns required to store the selected */
    /* eigenvectors, and standardize the array SELECT. */
    *m = 0;
    pair = FALSE_;
    i__1 = *n;
    for(k = 1; k <= i__1; ++k)
    {
        if(pair)
        {
            pair = FALSE_;
            select[k] = FALSE_;
        }
        else
        {
            if(wi[k] == 0.)
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
                    select[k] = TRUE_;
                    *m += 2;
                }
            }
        }
        /* L10: */
    }
    *info = 0;
    if(!rightv && !leftv)
    {
        *info = -1;
    }
    else if(!fromqr && !lsame_(eigsrc, "N", 1, 1))
    {
        *info = -2;
    }
    else if(!noinit && !lsame_(initv, "U", 1, 1))
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*ldh < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldvl < 1 || leftv && *ldvl < *n)
    {
        *info = -11;
    }
    else if(*ldvr < 1 || rightv && *ldvr < *n)
    {
        *info = -13;
    }
    else if(*mm < *m)
    {
        *info = -14;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DHSEIN", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible. */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Set machine-dependent constants. */
    unfl = dlamch_("Safe minimum");
    ulp = dlamch_("Precision");
    smlnum = unfl * (*n / ulp);
    bignum = (1. - ulp) / smlnum;
    eps3 = smlnum;
    ldwork = *n + 1;
    kl = 1;
    kln = 0;
    if(fromqr)
    {
        kr = 0;
    }
    else
    {
        kr = *n;
    }
    ksr = 1;
    i__1 = *n;
    for(k = 1; k <= i__1; ++k)
    {
        if(select[k])
        {
            /* Compute eigenvector(s) corresponding to W(K). */
            if(fromqr)
            {
                /* If affiliation of eigenvalues is known, check whether */
                /* the matrix splits. */
                /* Determine KL and KR such that 1 <= KL <= K <= KR <= N */
                /* and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or */
                /* KR = N). */
                /* Then inverse iteration can be performed with the */
                /* submatrix H(KL:N,KL:N) for a left eigenvector, and with */
                /* the submatrix H(1:KR,1:KR) for a right eigenvector. */
                i__2 = kl + 1;
                for(i__ = k; i__ >= i__2; --i__)
                {
                    if(h__[i__ + (i__ - 1) * h_dim1] == 0.)
                    {
                        goto L30;
                    }
                    /* L20: */
                }
            L30:
                kl = i__;
                if(k > kr)
                {
                    i__2 = *n - 1;
                    for(i__ = k; i__ <= i__2; ++i__)
                    {
                        if(h__[i__ + 1 + i__ * h_dim1] == 0.)
                        {
                            goto L50;
                        }
                        /* L40: */
                    }
                L50:
                    kr = i__;
                }
            }
            if(kl != kln)
            {
                kln = kl;
                /* Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it */
                /* has not ben computed before. */
                i__2 = kr - kl + 1;
                hnorm = aocl_lapack_dlanhs("I", &i__2, &h__[kl + kl * h_dim1], ldh, &work[1]);
                if(disnan_(&hnorm))
                {
                    *info = -6;
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
                else if(hnorm > 0.)
                {
                    eps3 = hnorm * ulp;
                }
                else
                {
                    eps3 = smlnum;
                }
            }
            /* Perturb eigenvalue if it is close to any previous */
            /* selected eigenvalues affiliated to the submatrix */
            /* H(KL:KR,KL:KR). Close roots are modified by EPS3. */
            wkr = wr[k];
            wki = wi[k];
        L60:
            i__2 = kl;
            for(i__ = k - 1; i__ >= i__2; --i__)
            {
                if(select[i__]
                   && (d__1 = wr[i__] - wkr, f2c_dabs(d__1))
                              + (d__2 = wi[i__] - wki, f2c_dabs(d__2))
                          < eps3)
                {
                    wkr += eps3;
                    goto L60;
                }
                /* L70: */
            }
            wr[k] = wkr;
            pair = wki != 0.;
            if(pair)
            {
                ksi = ksr + 1;
            }
            else
            {
                ksi = ksr;
            }
            if(leftv)
            {
                /* Compute left eigenvector. */
                i__2 = *n - kl + 1;
                aocl_lapack_dlaein(&c_false, &noinit, &i__2, &h__[kl + kl * h_dim1], ldh, &wkr,
                                   &wki, &vl[kl + ksr * vl_dim1], &vl[kl + ksi * vl_dim1], &work[1],
                                   &ldwork, &work[*n * *n + *n + 1], &eps3, &smlnum, &bignum,
                                   &iinfo);
                if(iinfo > 0)
                {
                    if(pair)
                    {
                        *info += 2;
                    }
                    else
                    {
                        ++(*info);
                    }
                    ifaill[ksr] = (aocl_int_t)(k);
                    ifaill[ksi] = (aocl_int_t)(k);
                }
                else
                {
                    ifaill[ksr] = 0;
                    ifaill[ksi] = 0;
                }
                i__2 = kl - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    vl[i__ + ksr * vl_dim1] = 0.;
                    /* L80: */
                }
                if(pair)
                {
                    i__2 = kl - 1;
                    for(i__ = 1; i__ <= i__2; ++i__)
                    {
                        vl[i__ + ksi * vl_dim1] = 0.;
                        /* L90: */
                    }
                }
            }
            if(rightv)
            {
                /* Compute right eigenvector. */
                aocl_lapack_dlaein(&c_true, &noinit, &kr, &h__[h_offset], ldh, &wkr, &wki,
                                   &vr[ksr * vr_dim1 + 1], &vr[ksi * vr_dim1 + 1], &work[1],
                                   &ldwork, &work[*n * *n + *n + 1], &eps3, &smlnum, &bignum,
                                   &iinfo);
                if(iinfo > 0)
                {
                    if(pair)
                    {
                        *info += 2;
                    }
                    else
                    {
                        ++(*info);
                    }
                    ifailr[ksr] = (aocl_int_t)(k);
                    ifailr[ksi] = (aocl_int_t)(k);
                }
                else
                {
                    ifailr[ksr] = 0;
                    ifailr[ksi] = 0;
                }
                i__2 = *n;
                for(i__ = kr + 1; i__ <= i__2; ++i__)
                {
                    vr[i__ + ksr * vr_dim1] = 0.;
                    /* L100: */
                }
                if(pair)
                {
                    i__2 = *n;
                    for(i__ = kr + 1; i__ <= i__2; ++i__)
                    {
                        vr[i__ + ksi * vr_dim1] = 0.;
                        /* L110: */
                    }
                }
            }
            if(pair)
            {
                ksr += 2;
            }
            else
            {
                ++ksr;
            }
        }
        /* L120: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DHSEIN */
}
/* dhsein_ */
