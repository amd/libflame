/* ./claqr2.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {0.f, 0.f};
static complex c_b2 = {1.f, 0.f};
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;
/* > \brief \b CLAQR2 performs the unitary similarity transformation of a Hessenberg matrix to
 * detect and defl ate fully converged eigenvalues from a trailing principal submatrix (aggressive
 * early deflation). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQR2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, */
/* IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, */
/* NV, WV, LDWV, WORK, LWORK ) */
/* .. Scalar Arguments .. */
/* INTEGER IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, */
/* $ LDZ, LWORK, N, ND, NH, NS, NV, NW */
/* LOGICAL WANTT, WANTZ */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), */
/* $ WORK( * ), WV( LDWV, * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQR2 is identical to CLAQR3 except that it avoids */
/* > recursion by calling CLAHQR instead of CLAQR4. */
/* > */
/* > Aggressive early deflation: */
/* > */
/* > This subroutine accepts as input an upper Hessenberg matrix */
/* > H and performs an unitary similarity transformation */
/* > designed to detect and deflate fully converged eigenvalues from */
/* > a trailing principal submatrix. On output H has been over- */
/* > written by a new Hessenberg matrix that is a perturbation of */
/* > an unitary similarity transformation of H. It is to be */
/* > hoped that the final version of H has many zero subdiagonal */
/* > entries. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTT */
/* > \verbatim */
/* > WANTT is LOGICAL */
/* > If .TRUE., then the Hessenberg matrix H is fully updated */
/* > so that the triangular Schur factor may be */
/* > computed (in cooperation with the calling subroutine). */
/* > If .FALSE., then only enough of H is updated to preserve */
/* > the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > If .TRUE., then the unitary matrix Z is updated so */
/* > so that the unitary Schur factor may be computed */
/* > (in cooperation with the calling subroutine). */
/* > If .FALSE., then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix H and (if WANTZ is .TRUE.) the */
/* > order of the unitary matrix Z. */
/* > \endverbatim */
/* > */
/* > \param[in] KTOP */
/* > \verbatim */
/* > KTOP is INTEGER */
/* > It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0. */
/* > KBOT and KTOP together determine an isolated block */
/* > along the diagonal of the Hessenberg matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] KBOT */
/* > \verbatim */
/* > KBOT is INTEGER */
/* > It is assumed without a check that either */
/* > KBOT = N or H(KBOT+1,KBOT)=0. KBOT and KTOP together */
/* > determine an isolated block along the diagonal of the */
/* > Hessenberg matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] NW */
/* > \verbatim */
/* > NW is INTEGER */
/* > Deflation window size. 1 <= NW <= (KBOT-KTOP+1). */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* > H is COMPLEX array, dimension (LDH,N) */
/* > On input the initial N-by-N section of H stores the */
/* > Hessenberg matrix undergoing aggressive early deflation. */
/* > On output H has been transformed by a unitary */
/* > similarity transformation, perturbed, and the returned */
/* > to Hessenberg form that (it is to be hoped) has some */
/* > zero subdiagonal entries. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > Leading dimension of H just as declared in the calling */
/* > subroutine. N <= LDH */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* > ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* > IHIZ is INTEGER */
/* > Specify the rows of Z to which transformations must be */
/* > applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ,N) */
/* > IF WANTZ is .TRUE., then on output, the unitary */
/* > similarity transformation mentioned above has been */
/* > accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right. */
/* > If WANTZ is .FALSE., then Z is unreferenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of Z just as declared in the */
/* > calling subroutine. 1 <= LDZ. */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* > NS is INTEGER */
/* > The number of unconverged (ie approximate) eigenvalues */
/* > returned in SR and SI that may be used as shifts by the */
/* > calling subroutine. */
/* > \endverbatim */
/* > */
/* > \param[out] ND */
/* > \verbatim */
/* > ND is INTEGER */
/* > The number of converged eigenvalues uncovered by this */
/* > subroutine. */
/* > \endverbatim */
/* > */
/* > \param[out] SH */
/* > \verbatim */
/* > SH is COMPLEX array, dimension (KBOT) */
/* > On output, approximate eigenvalues that may */
/* > be used for shifts are stored in SH(KBOT-ND-NS+1) */
/* > through SR(KBOT-ND). Converged eigenvalues are */
/* > stored in SH(KBOT-ND+1) through SH(KBOT). */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension (LDV,NW) */
/* > An NW-by-NW work array. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of V just as declared in the */
/* > calling subroutine. NW <= LDV */
/* > \endverbatim */
/* > */
/* > \param[in] NH */
/* > \verbatim */
/* > NH is INTEGER */
/* > The number of columns of T. NH >= NW. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX array, dimension (LDT,NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of T just as declared in the */
/* > calling subroutine. NW <= LDT */
/* > \endverbatim */
/* > */
/* > \param[in] NV */
/* > \verbatim */
/* > NV is INTEGER */
/* > The number of rows of work array WV available for */
/* > workspace. NV >= NW. */
/* > \endverbatim */
/* > */
/* > \param[out] WV */
/* > \verbatim */
/* > WV is COMPLEX array, dimension (LDWV,NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWV */
/* > \verbatim */
/* > LDWV is INTEGER */
/* > The leading dimension of W just as declared in the */
/* > calling subroutine. NW <= LDV */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (LWORK) */
/* > On exit, WORK(1) is set to an estimate of the optimal value */
/* > of LWORK for the given values of N, NW, KTOP and KBOT. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the work array WORK. LWORK = 2*NW */
/* > suffices, but greater efficiency may result from larger */
/* > values of LWORK. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
CLAQR2 */
/* > only estimates the optimal workspace size for the given */
/* > values of N, NW, KTOP and KBOT. The estimate is returned */
/* > in WORK(1). No error message related to LWORK is issued */
/* > by XERBLA. Neither H nor Z are accessed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup laqr2 */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Karen Braman and Ralph Byers, Department of Mathematics, */
/* > University of Kansas, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
void claqr2_(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw,
             complex *h__, integer *ldh, integer *iloz, integer *ihiz, complex *z__, integer *ldz,
             integer *ns, integer *nd, complex *sh, complex *v, integer *ldv, integer *nh,
             complex *t, integer *ldt, integer *nv, complex *wv, integer *ldwv, complex *work,
             integer *lwork)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "claqr2 inputs: n %lld, ktop %lld, kbot %lld, nw %lld, ldh %lld, iloz %lld, ihiz "
             "%lld, ldz %lld, ldv %lld, nh %lld, ldt %lld, nv %lld, ldwv %lld, lwork %lld",
             *n, *ktop, *kbot, *nw, *ldh, *iloz, *ihiz, *ldz, *ldv, *nh, *ldt, *nv, *ldwv, *lwork);
#else
    snprintf(buffer, 256,
             "claqr2 inputs: n %d, ktop %d, kbot %d, nw %d, ldh %d, iloz %d, ihiz %d, ldz %d, ldv "
             "%d, nh %d, ldt %d, nv %d, ldwv %d, lwork %d",
             *n, *ktop, *kbot, *nw, *ldh, *iloz, *ihiz, *ldz, *ldv, *nh, *ldt, *nv, *ldwv, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer h_dim1, h_offset, t_dim1, t_offset, v_dim1, v_offset, wv_dim1, wv_offset, z_dim1,
        z_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer i__, j;
    complex s;
    integer jw;
    real foo;
    integer kln;
    complex tau;
    integer knt;
    real ulp;
    integer lwk1, lwk2;
    complex beta;
    integer kcol, info, ifst, ilst, ltop, krow;
    extern /* Subroutine */
        void
        clarf_(char *, integer *, integer *, complex *, integer *, complex *, complex *, integer *,
               complex *),
        cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, complex *, integer *),
        ccopy_(integer *, complex *, integer *, complex *, integer *);
    integer infqr, kwtop;
    extern /* Subroutine */
        void
        cgehrd_(integer *, integer *, integer *, complex *, integer *, complex *, complex *,
                integer *, integer *),
        clarfg_(integer *, complex *, complex *, integer *, complex *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        clahqr_(logical *, logical *, integer *, integer *, integer *, complex *, integer *,
                complex *, integer *, integer *, complex *, integer *, integer *),
        clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *),
        claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *);
    real safmin;
    extern /* Subroutine */
        void
        ctrexc_(char *, integer *, complex *, integer *, complex *, integer *, integer *, integer *,
                integer *),
        cunmhr_(char *, char *, integer *, integer *, integer *, integer *, complex *, integer *,
                complex *, complex *, integer *, complex *, integer *, integer *);
    real smlnum;
    integer lwkopt;
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ================================================================ */
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* ==== Estimate optimal workspace. ==== */
    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --sh;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    wv_dim1 = *ldwv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    --work;
    /* Function Body */
    /* Computing MIN */
    i__1 = *nw;
    i__2 = *kbot - *ktop + 1; // , expr subst
    jw = fla_min(i__1, i__2);
    if(jw <= 2)
    {
        lwkopt = 1;
    }
    else
    {
        /* ==== Workspace query call to CGEHRD ==== */
        i__1 = jw - 1;
        cgehrd_(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &c_n1, &info);
        lwk1 = (integer)work[1].r;
        /* ==== Workspace query call to CUNMHR ==== */
        i__1 = jw - 1;
        cunmhr_("R", "N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &v[v_offset], ldv,
                &work[1], &c_n1, &info);
        lwk2 = (integer)work[1].r;
        /* ==== Optimal workspace ==== */
        lwkopt = jw + fla_max(lwk1, lwk2);
    }
    /* ==== Quick return in case of workspace query. ==== */
    if(*lwork == -1)
    {
        r__1 = (real)lwkopt;
        q__1.r = r__1;
        q__1.i = 0.f; // , expr subst
        work[1].r = q__1.r;
        work[1].i = q__1.i; // , expr subst
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* ==== Nothing to do ... */
    /* ... for an empty active block ... ==== */
    *ns = 0;
    *nd = 0;
    work[1].r = 1.f;
    work[1].i = 0.f; // , expr subst
    if(*ktop > *kbot)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* ... nor for an empty deflation window. ==== */
    if(*nw < 1)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* ==== Machine constants ==== */
    safmin = slamch_("SAFE MINIMUM");
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real)(*n) / ulp);
    /* ==== Setup deflation window ==== */
    /* Computing MIN */
    i__1 = *nw;
    i__2 = *kbot - *ktop + 1; // , expr subst
    jw = fla_min(i__1, i__2);
    kwtop = *kbot - jw + 1;
    if(kwtop == *ktop)
    {
        s.r = 0.f;
        s.i = 0.f; // , expr subst
    }
    else
    {
        i__1 = kwtop + (kwtop - 1) * h_dim1;
        s.r = h__[i__1].r;
        s.i = h__[i__1].i; // , expr subst
    }
    if(*kbot == kwtop)
    {
        /* ==== 1-by-1 deflation window: not much to do ==== */
        i__1 = kwtop;
        i__2 = kwtop + kwtop * h_dim1;
        sh[i__1].r = h__[i__2].r;
        sh[i__1].i = h__[i__2].i; // , expr subst
        *ns = 1;
        *nd = 0;
        /* Computing MAX */
        i__1 = kwtop + kwtop * h_dim1;
        r__5 = smlnum;
        r__6 = ulp
               * ((r__1 = h__[i__1].r, f2c_abs(r__1))
                  + (r__2 = r_imag(&h__[kwtop + kwtop * h_dim1]), f2c_abs(r__2))); // , expr subst
        if((r__3 = s.r, f2c_abs(r__3)) + (r__4 = r_imag(&s), f2c_abs(r__4)) <= fla_max(r__5, r__6))
        {
            *ns = 0;
            *nd = 1;
            if(kwtop > *ktop)
            {
                i__1 = kwtop + (kwtop - 1) * h_dim1;
                h__[i__1].r = 0.f;
                h__[i__1].i = 0.f; // , expr subst
            }
        }
        work[1].r = 1.f;
        work[1].i = 0.f; // , expr subst
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* ==== Convert to spike-triangular form. (In case of a */
    /* . rare QR failure, this routine continues to do */
    /* . aggressive early deflation using that part of */
    /* . the deflation window that converged using INFQR */
    /* . here and there to keep track.) ==== */
    clacpy_("U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], ldt);
    i__1 = jw - 1;
    i__2 = *ldh + 1;
    i__3 = *ldt + 1;
    ccopy_(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &i__3);
    claset_("A", &jw, &jw, &c_b1, &c_b2, &v[v_offset], ldv);
    clahqr_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sh[kwtop], &c__1, &jw,
            &v[v_offset], ldv, &infqr);
    /* ==== Deflation detection loop ==== */
    *ns = jw;
    ilst = infqr + 1;
    i__1 = jw;
    for(knt = infqr + 1; knt <= i__1; ++knt)
    {
        /* ==== Small spike tip deflation test ==== */
        i__2 = *ns + *ns * t_dim1;
        foo = (r__1 = t[i__2].r, f2c_abs(r__1))
              + (r__2 = r_imag(&t[*ns + *ns * t_dim1]), f2c_abs(r__2));
        if(foo == 0.f)
        {
            foo = (r__1 = s.r, f2c_abs(r__1)) + (r__2 = r_imag(&s), f2c_abs(r__2));
        }
        i__2 = *ns * v_dim1 + 1;
        /* Computing MAX */
        r__5 = smlnum;
        r__6 = ulp * foo; // , expr subst
        if(((r__1 = s.r, f2c_abs(r__1)) + (r__2 = r_imag(&s), f2c_abs(r__2)))
               * ((r__3 = v[i__2].r, f2c_abs(r__3))
                  + (r__4 = r_imag(&v[*ns * v_dim1 + 1]), f2c_abs(r__4)))
           <= fla_max(r__5, r__6))
        {
            /* ==== One more converged eigenvalue ==== */
            --(*ns);
        }
        else
        {
            /* ==== One undeflatable eigenvalue. Move it up out of the */
            /* . way. (CTREXC can not fail in this case.) ==== */
            ifst = *ns;
            ctrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &ilst, &info);
            ++ilst;
        }
        /* L10: */
    }
    /* ==== Return to Hessenberg form ==== */
    if(*ns == 0)
    {
        s.r = 0.f;
        s.i = 0.f; // , expr subst
    }
    if(*ns < jw)
    {
        /* ==== sorting the diagonal of T improves accuracy for */
        /* . graded matrices. ==== */
        i__1 = *ns;
        for(i__ = infqr + 1; i__ <= i__1; ++i__)
        {
            ifst = i__;
            i__2 = *ns;
            for(j = i__ + 1; j <= i__2; ++j)
            {
                i__3 = j + j * t_dim1;
                i__4 = ifst + ifst * t_dim1;
                if((r__1 = t[i__3].r, f2c_abs(r__1))
                       + (r__2 = r_imag(&t[j + j * t_dim1]), f2c_abs(r__2))
                   > (r__3 = t[i__4].r, f2c_abs(r__3))
                         + (r__4 = r_imag(&t[ifst + ifst * t_dim1]), f2c_abs(r__4)))
                {
                    ifst = j;
                }
                /* L20: */
            }
            ilst = i__;
            if(ifst != ilst)
            {
                ctrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &ilst, &info);
            }
            /* L30: */
        }
    }
    /* ==== Restore shift/eigenvalue array from T ==== */
    i__1 = jw;
    for(i__ = infqr + 1; i__ <= i__1; ++i__)
    {
        i__2 = kwtop + i__ - 1;
        i__3 = i__ + i__ * t_dim1;
        sh[i__2].r = t[i__3].r;
        sh[i__2].i = t[i__3].i; // , expr subst
        /* L40: */
    }
    if(*ns < jw || s.r == 0.f && s.i == 0.f)
    {
        if(*ns > 1 && (s.r != 0.f || s.i != 0.f))
        {
            /* ==== Reflect spike back into lower triangle ==== */
            ccopy_(ns, &v[v_offset], ldv, &work[1], &c__1);
            i__1 = *ns;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                r_cnjg(&q__1, &work[i__]);
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
                /* L50: */
            }
            beta.r = work[1].r;
            beta.i = work[1].i; // , expr subst
            clarfg_(ns, &beta, &work[2], &c__1, &tau);
            work[1].r = 1.f;
            work[1].i = 0.f; // , expr subst
            i__1 = jw - 2;
            i__2 = jw - 2;
            claset_("L", &i__1, &i__2, &c_b1, &c_b1, &t[t_dim1 + 3], ldt);
            r_cnjg(&q__1, &tau);
            clarf_("L", ns, &jw, &work[1], &c__1, &q__1, &t[t_offset], ldt, &work[jw + 1]);
            clarf_("R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &work[jw + 1]);
            clarf_("R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv, &work[jw + 1]);
            i__1 = *lwork - jw;
            cgehrd_(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1], &i__1, &info);
        }
        /* ==== Copy updated reduced window into place ==== */
        if(kwtop > 1)
        {
            i__1 = kwtop + (kwtop - 1) * h_dim1;
            r_cnjg(&q__2, &v[v_dim1 + 1]);
            q__1.r = s.r * q__2.r - s.i * q__2.i;
            q__1.i = s.r * q__2.i + s.i * q__2.r; // , expr subst
            h__[i__1].r = q__1.r;
            h__[i__1].i = q__1.i; // , expr subst
        }
        clacpy_("U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1], ldh);
        i__1 = jw - 1;
        i__2 = *ldt + 1;
        i__3 = *ldh + 1;
        ccopy_(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1], &i__3);
        /* ==== Accumulate orthogonal matrix in order update */
        /* . H and Z, if requested. ==== */
        if(*ns > 1 && (s.r != 0.f || s.i != 0.f))
        {
            i__1 = *lwork - jw;
            cunmhr_("R", "N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1], &v[v_offset], ldv,
                    &work[jw + 1], &i__1, &info);
        }
        /* ==== Update vertical slab in H ==== */
        if(*wantt)
        {
            ltop = 1;
        }
        else
        {
            ltop = *ktop;
        }
        i__1 = kwtop - 1;
        i__2 = *nv;
        for(krow = ltop; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += i__2)
        {
            /* Computing MIN */
            i__3 = *nv;
            i__4 = kwtop - krow; // , expr subst
            kln = fla_min(i__3, i__4);
            cgemm_("N", "N", &kln, &jw, &jw, &c_b2, &h__[krow + kwtop * h_dim1], ldh, &v[v_offset],
                   ldv, &c_b1, &wv[wv_offset], ldwv);
            clacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * h_dim1], ldh);
            /* L60: */
        }
        /* ==== Update horizontal slab in H ==== */
        if(*wantt)
        {
            i__2 = *n;
            i__1 = *nh;
            for(kcol = *kbot + 1; i__1 < 0 ? kcol >= i__2 : kcol <= i__2; kcol += i__1)
            {
                /* Computing MIN */
                i__3 = *nh;
                i__4 = *n - kcol + 1; // , expr subst
                kln = fla_min(i__3, i__4);
                cgemm_("C", "N", &jw, &kln, &jw, &c_b2, &v[v_offset], ldv,
                       &h__[kwtop + kcol * h_dim1], ldh, &c_b1, &t[t_offset], ldt);
                clacpy_("A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol * h_dim1], ldh);
                /* L70: */
            }
        }
        /* ==== Update vertical slab in Z ==== */
        if(*wantz)
        {
            i__1 = *ihiz;
            i__2 = *nv;
            for(krow = *iloz; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += i__2)
            {
                /* Computing MIN */
                i__3 = *nv;
                i__4 = *ihiz - krow + 1; // , expr subst
                kln = fla_min(i__3, i__4);
                cgemm_("N", "N", &kln, &jw, &jw, &c_b2, &z__[krow + kwtop * z_dim1], ldz,
                       &v[v_offset], ldv, &c_b1, &wv[wv_offset], ldwv);
                clacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &z__[krow + kwtop * z_dim1], ldz);
                /* L80: */
            }
        }
    }
    /* ==== Return the number of deflations ... ==== */
    *nd = jw - *ns;
    /* ==== ... and the number of shifts. (Subtracting */
    /* . INFQR from the spike length takes care */
    /* . of the case of a rare QR failure while */
    /* . calculating eigenvalues of the deflation */
    /* . window.) ==== */
    *ns -= infqr;
    /* ==== Return optimal workspace. ==== */
    r__1 = (real)lwkopt;
    q__1.r = r__1;
    q__1.i = 0.f; // , expr subst
    work[1].r = q__1.r;
    work[1].i = q__1.i; // , expr subst
    /* ==== End of CLAQR2 ==== */
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* claqr2_ */
