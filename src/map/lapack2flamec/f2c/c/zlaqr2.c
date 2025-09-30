/* ./zlaqr2.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {{0.}, {0.}};
static dcomplex c_b2 = {{1.}, {0.}};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static logical c_true = TRUE_;
/* > \brief \b ZLAQR2 performs the unitary similarity transformation of a Hessenberg matrix to
 * detect and defl ate fully converged eigenvalues from a trailing principal submatrix (aggressive
 * early deflation). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAQR2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, */
/* IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, */
/* NV, WV, LDWV, WORK, LWORK ) */
/* .. Scalar Arguments .. */
/* INTEGER IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, */
/* $ LDZ, LWORK, N, ND, NH, NS, NV, NW */
/* LOGICAL WANTT, WANTZ */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), */
/* $ WORK( * ), WV( LDWV, * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQR2 is identical to ZLAQR3 except that it avoids */
/* > recursion by calling ZLAHQR instead of ZLAQR4. */
/* > */
/* > Aggressive early deflation: */
/* > */
/* > ZLAQR2 accepts as input an upper Hessenberg matrix */
/* > H and performs an unitary similarity transformation */
/* > designed to detect and deflate fully converged eigenvalues from */
/* > a trailing principal submatrix. On output H has been over- */
/* > written by a new Hessenberg matrix that is a perturbation of */
/* > an unitary similarity transformation of H. It is to be */
/* > hoped that the final version of H has many zero subdiagonal */
/* > entries. */
/* > */
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
/* > H is COMPLEX*16 array, dimension (LDH,N) */
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
/* > Z is COMPLEX*16 array, dimension (LDZ,N) */
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
/* > SH is COMPLEX*16 array, dimension (KBOT) */
/* > On output, approximate eigenvalues that may */
/* > be used for shifts are stored in SH(KBOT-ND-NS+1) */
/* > through SR(KBOT-ND). Converged eigenvalues are */
/* > stored in SH(KBOT-ND+1) through SH(KBOT). */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX*16 array, dimension (LDV,NW) */
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
/* > T is COMPLEX*16 array, dimension (LDT,NW) */
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
/* > WV is COMPLEX*16 array, dimension (LDWV,NW) */
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
/* > WORK is COMPLEX*16 array, dimension (LWORK) */
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
ZLAQR2 */
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
/** Generated wrapper function */
void zlaqr2_(logical *wantt, logical *wantz, aocl_int_t *n, aocl_int_t *ktop, aocl_int_t *kbot,
             aocl_int_t *nw, dcomplex *h__, aocl_int_t *ldh, aocl_int_t *iloz,
             aocl_int_t *ihiz, dcomplex *z__, aocl_int_t *ldz, aocl_int_t *ns, aocl_int_t *nd,
             dcomplex *sh, dcomplex *v, aocl_int_t *ldv, aocl_int_t *nh, dcomplex *t,
             aocl_int_t *ldt, aocl_int_t *nv, dcomplex *wv, aocl_int_t *ldwv,
             dcomplex *work, aocl_int_t *lwork)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlaqr2(wantt, wantz, n, ktop, kbot, nw, h__, ldh, iloz, ihiz, z__, ldz, ns, nd, sh,
                       v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ktop_64 = *ktop;
    aocl_int64_t kbot_64 = *kbot;
    aocl_int64_t nw_64 = *nw;
    aocl_int64_t ldh_64 = *ldh;
    aocl_int64_t iloz_64 = *iloz;
    aocl_int64_t ihiz_64 = *ihiz;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t ns_64 = *ns;
    aocl_int64_t nd_64 = *nd;
    aocl_int64_t ldv_64 = *ldv;
    aocl_int64_t nh_64 = *nh;
    aocl_int64_t ldt_64 = *ldt;
    aocl_int64_t nv_64 = *nv;
    aocl_int64_t ldwv_64 = *ldwv;
    aocl_int64_t lwork_64 = *lwork;

    aocl_lapack_zlaqr2(wantt, wantz, &n_64, &ktop_64, &kbot_64, &nw_64, h__, &ldh_64, &iloz_64,
                       &ihiz_64, z__, &ldz_64, &ns_64, &nd_64, sh, v, &ldv_64, &nh_64, t, &ldt_64,
                       &nv_64, wv, &ldwv_64, work, &lwork_64);

    *ns = (aocl_int_t)ns_64;
    *nd = (aocl_int_t)nd_64;
#endif
}

void aocl_lapack_zlaqr2(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop,
                        aocl_int64_t *kbot, aocl_int64_t *nw, dcomplex *h__, aocl_int64_t *ldh,
                        aocl_int64_t *iloz, aocl_int64_t *ihiz, dcomplex *z__,
                        aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, dcomplex *sh,
                        dcomplex *v, aocl_int64_t *ldv, aocl_int64_t *nh, dcomplex *t,
                        aocl_int64_t *ldt, aocl_int64_t *nv, dcomplex *wv, aocl_int64_t *ldwv,
                        dcomplex *work, aocl_int64_t *lwork)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
        "zlaqr2 inputs: n %" FLA_IS ", ktop %" FLA_IS ", kbot %" FLA_IS ", nw %" FLA_IS
        ", ldh %" FLA_IS ", iloz %" FLA_IS ", ihiz %" FLA_IS ", ldz %" FLA_IS ", ldv %" FLA_IS
        ", nh %" FLA_IS ", ldt %" FLA_IS ", nv %" FLA_IS ", ldwv %" FLA_IS ",lwork %" FLA_IS "",
        *n, *ktop, *kbot, *nw, *ldh, *iloz, *ihiz, *ldz, *ldv, *nh, *ldt, *nv, *ldwv, *lwork);
    /* System generated locals */
    aocl_int64_t h_dim1, h_offset, t_dim1, t_offset, v_dim1, v_offset, wv_dim1, wv_offset, z_dim1,
        z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    dcomplex z__1, z__2;
    /* Builtin functions */
    double d_imag(dcomplex *);
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__, j;
    dcomplex s;
    aocl_int64_t jw;
    doublereal foo;
    aocl_int64_t kln;
    dcomplex tau;
    aocl_int64_t knt;
    doublereal ulp;
    aocl_int64_t lwk1, lwk2;
    dcomplex beta;
    aocl_int64_t kcol, info, ifst, ilst, ltop, krow;
    aocl_int64_t infqr;
    aocl_int64_t kwtop;
    extern doublereal dlamch_(char *);
    doublereal safmin;
    doublereal smlnum;
    aocl_int64_t lwkopt;
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
        /* ==== Workspace query call to ZGEHRD ==== */
        i__1 = jw - 1;
        aocl_lapack_zgehrd(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &c_n1, &info);
        lwk1 = (integer)work[1].r;
        /* ==== Workspace query call to ZUNMHR ==== */
        i__1 = jw - 1;
        aocl_lapack_zunmhr("R", "N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1],
                           &v[v_offset], ldv, &work[1], &c_n1, &info);
        lwk2 = (integer)work[1].r;
        /* ==== Optimal workspace ==== */
        lwkopt = jw + fla_max(lwk1, lwk2);
    }
    /* ==== Quick return in case of workspace query. ==== */
    if(*lwork == -1)
    {
        d__1 = (doublereal)lwkopt;
        z__1.r = d__1;
        z__1.i = 0.; // , expr subst
        work[1].r = z__1.r;
        work[1].i = z__1.i; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ==== Nothing to do ... */
    /* ... for an empty active block ... ==== */
    *ns = 0;
    *nd = 0;
    work[1].r = 1.;
    work[1].i = 0.; // , expr subst
    if(*ktop > *kbot)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ... nor for an empty deflation window. ==== */
    if(*nw < 1)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ==== Machine constants ==== */
    safmin = dlamch_("SAFE MINIMUM");
    ulp = dlamch_("PRECISION");
    smlnum = safmin * ((doublereal)(*n) / ulp);
    /* ==== Setup deflation window ==== */
    /* Computing MIN */
    i__1 = *nw;
    i__2 = *kbot - *ktop + 1; // , expr subst
    jw = fla_min(i__1, i__2);
    kwtop = *kbot - jw + 1;
    if(kwtop == *ktop)
    {
        s.r = 0.;
        s.i = 0.; // , expr subst
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
        d__5 = smlnum;
        d__6 = ulp
               * ((d__1 = h__[i__1].r, f2c_dabs(d__1))
                  + (d__2 = d_imag(&h__[kwtop + kwtop * h_dim1]), f2c_dabs(d__2))); // , expr subst
        if((d__3 = s.r, f2c_dabs(d__3)) + (d__4 = d_imag(&s), f2c_dabs(d__4))
           <= fla_max(d__5, d__6))
        {
            *ns = 0;
            *nd = 1;
            if(kwtop > *ktop)
            {
                i__1 = kwtop + (kwtop - 1) * h_dim1;
                h__[i__1].r = 0.;
                h__[i__1].i = 0.; // , expr subst
            }
        }
        work[1].r = 1.;
        work[1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ==== Convert to spike-triangular form. (In case of a */
    /* . rare QR failure, this routine continues to do */
    /* . aggressive early deflation using that part of */
    /* . the deflation window that converged using INFQR */
    /* . here and there to keep track.) ==== */
    aocl_lapack_zlacpy("U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], ldt);
    i__1 = jw - 1;
    i__2 = *ldh + 1;
    i__3 = *ldt + 1;
    aocl_blas_zcopy(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &i__3);
    aocl_lapack_zlaset("A", &jw, &jw, &c_b1, &c_b2, &v[v_offset], ldv);
    aocl_lapack_zlahqr(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sh[kwtop], &c__1, &jw,
                       &v[v_offset], ldv, &infqr);
    /* ==== Deflation detection loop ==== */
    *ns = jw;
    ilst = infqr + 1;
    i__1 = jw;
    for(knt = infqr + 1; knt <= i__1; ++knt)
    {
        /* ==== Small spike tip deflation test ==== */
        i__2 = *ns + *ns * t_dim1;
        foo = (d__1 = t[i__2].r, f2c_dabs(d__1))
              + (d__2 = d_imag(&t[*ns + *ns * t_dim1]), f2c_dabs(d__2));
        if(foo == 0.)
        {
            foo = (d__1 = s.r, f2c_dabs(d__1)) + (d__2 = d_imag(&s), f2c_dabs(d__2));
        }
        i__2 = *ns * v_dim1 + 1;
        /* Computing MAX */
        d__5 = smlnum;
        d__6 = ulp * foo; // , expr subst
        if(((d__1 = s.r, f2c_dabs(d__1)) + (d__2 = d_imag(&s), f2c_dabs(d__2)))
               * ((d__3 = v[i__2].r, f2c_dabs(d__3))
                  + (d__4 = d_imag(&v[*ns * v_dim1 + 1]), f2c_dabs(d__4)))
           <= fla_max(d__5, d__6))
        {
            /* ==== One more converged eigenvalue ==== */
            --(*ns);
        }
        else
        {
            /* ==== One undeflatable eigenvalue. Move it up out of the */
            /* . way. (ZTREXC can not fail in this case.) ==== */
            ifst = *ns;
            aocl_lapack_ztrexc("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &ilst, &info);
            ++ilst;
        }
        /* L10: */
    }
    /* ==== Return to Hessenberg form ==== */
    if(*ns == 0)
    {
        s.r = 0.;
        s.i = 0.; // , expr subst
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
                if((d__1 = t[i__3].r, f2c_dabs(d__1))
                       + (d__2 = d_imag(&t[j + j * t_dim1]), f2c_dabs(d__2))
                   > (d__3 = t[i__4].r, f2c_dabs(d__3))
                         + (d__4 = d_imag(&t[ifst + ifst * t_dim1]), f2c_dabs(d__4)))
                {
                    ifst = j;
                }
                /* L20: */
            }
            ilst = i__;
            if(ifst != ilst)
            {
                aocl_lapack_ztrexc("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &ilst,
                                   &info);
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
    if(*ns < jw || s.r == 0. && s.i == 0.)
    {
        if(*ns > 1 && (s.r != 0. || s.i != 0.))
        {
            /* ==== Reflect spike back into lower triangle ==== */
            aocl_blas_zcopy(ns, &v[v_offset], ldv, &work[1], &c__1);
            i__1 = *ns;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                d_cnjg(&z__1, &work[i__]);
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
                /* L50: */
            }
            beta.r = work[1].r;
            beta.i = work[1].i; // , expr subst
            aocl_lapack_zlarfg(ns, &beta, &work[2], &c__1, &tau);
            work[1].r = 1.;
            work[1].i = 0.; // , expr subst
            i__1 = jw - 2;
            i__2 = jw - 2;
            aocl_lapack_zlaset("L", &i__1, &i__2, &c_b1, &c_b1, &t[t_dim1 + 3], ldt);
            d_cnjg(&z__1, &tau);
            aocl_lapack_zlarf("L", ns, &jw, &work[1], &c__1, &z__1, &t[t_offset], ldt,
                              &work[jw + 1]);
            aocl_lapack_zlarf("R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &work[jw + 1]);
            aocl_lapack_zlarf("R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv,
                              &work[jw + 1]);
            i__1 = *lwork - jw;
            aocl_lapack_zgehrd(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1], &i__1,
                               &info);
        }
        /* ==== Copy updated reduced window into place ==== */
        if(kwtop > 1)
        {
            i__1 = kwtop + (kwtop - 1) * h_dim1;
            d_cnjg(&z__2, &v[v_dim1 + 1]);
            z__1.r = s.r * z__2.r - s.i * z__2.i;
            z__1.i = s.r * z__2.i + s.i * z__2.r; // , expr subst
            h__[i__1].r = z__1.r;
            h__[i__1].i = z__1.i; // , expr subst
        }
        aocl_lapack_zlacpy("U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1], ldh);
        i__1 = jw - 1;
        i__2 = *ldt + 1;
        i__3 = *ldh + 1;
        aocl_blas_zcopy(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1], &i__3);
        /* ==== Accumulate orthogonal matrix in order update */
        /* . H and Z, if requested. ==== */
        if(*ns > 1 && (s.r != 0. || s.i != 0.))
        {
            i__1 = *lwork - jw;
            aocl_lapack_zunmhr("R", "N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1],
                               &v[v_offset], ldv, &work[jw + 1], &i__1, &info);
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
            aocl_blas_zgemm("N", "N", &kln, &jw, &jw, &c_b2, &h__[krow + kwtop * h_dim1], ldh,
                            &v[v_offset], ldv, &c_b1, &wv[wv_offset], ldwv);
            aocl_lapack_zlacpy("A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * h_dim1],
                               ldh);
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
                aocl_blas_zgemm("C", "N", &jw, &kln, &jw, &c_b2, &v[v_offset], ldv,
                                &h__[kwtop + kcol * h_dim1], ldh, &c_b1, &t[t_offset], ldt);
                aocl_lapack_zlacpy("A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol * h_dim1],
                                   ldh);
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
                aocl_blas_zgemm("N", "N", &kln, &jw, &jw, &c_b2, &z__[krow + kwtop * z_dim1], ldz,
                                &v[v_offset], ldv, &c_b1, &wv[wv_offset], ldwv);
                aocl_lapack_zlacpy("A", &kln, &jw, &wv[wv_offset], ldwv,
                                   &z__[krow + kwtop * z_dim1], ldz);
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
    d__1 = (doublereal)lwkopt;
    z__1.r = d__1;
    z__1.i = 0.; // , expr subst
    work[1].r = z__1.r;
    work[1].i = z__1.i; // , expr subst
    /* ==== End of ZLAQR2 ==== */
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* zlaqr2_ */
