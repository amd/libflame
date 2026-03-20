/* claqr5.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {0.f, 0.f};
static scomplex c_b2 = {1.f, 0.f};
static aocl_int64_t c__2 = 2;
static aocl_int64_t c__1 = 1;
static aocl_int64_t c__3 = 3;
/* > \brief \b CLAQR5 performs a single small-bulge multi-shift QR sweep. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQR5 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr5.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr5.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr5.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S, */
/* H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, */
/* WV, LDWV, NH, WH, LDWH ) */
/* .. Scalar Arguments .. */
/* INTEGER IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, */
/* $ LDWH, LDWV, LDZ, N, NH, NSHFTS, NV */
/* LOGICAL WANTT, WANTZ */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ), */
/* $ WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQR5 called by CLAQR0 performs a */
/* > single small-bulge multi-shift QR sweep. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTT */
/* > \verbatim */
/* > WANTT is LOGICAL */
/* > WANTT = .true. if the triangular Schur factor */
/* > is being computed. WANTT is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > WANTZ = .true. if the unitary Schur factor is being */
/* > computed. WANTZ is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] KACC22 */
/* > \verbatim */
/* > KACC22 is INTEGER with value 0, 1, or 2. */
/* > Specifies the computation mode of far-from-diagonal */
/* > orthogonal updates. */
/* > = 0: CLAQR5 does not accumulate reflections and does not */
/* > use matrix-matrix multiply to update far-from-diagonal */
/* > matrix entries. */
/* > = 1: CLAQR5 accumulates reflections and uses matrix-matrix */
/* > multiply to update the far-from-diagonal matrix entries. */
/* > = 2: Same as KACC22 = 1. This option used to enable exploiting */
/* > the 2-by-2 structure during matrix multiplications, but */
/* > this is no longer supported. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > N is the order of the Hessenberg matrix H upon which this */
/* > subroutine operates. */
/* > \endverbatim */
/* > */
/* > \param[in] KTOP */
/* > \verbatim */
/* > KTOP is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] KBOT */
/* > \verbatim */
/* > KBOT is INTEGER */
/* > These are the first and last rows and columns of an */
/* > isolated diagonal block upon which the QR sweep is to be */
/* > applied. It is assumed without a check that */
/* > either KTOP = 1 or H(KTOP,KTOP-1) = 0 */
/* > and */
/* > either KBOT = N or H(KBOT+1,KBOT) = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NSHFTS */
/* > \verbatim */
/* > NSHFTS is INTEGER */
/* > NSHFTS gives the number of simultaneous shifts. NSHFTS */
/* > must be positive and even. */
/* > \endverbatim */
/* > */
/* > \param[in,out] S */
/* > \verbatim */
/* > S is COMPLEX array, dimension (NSHFTS) */
/* > S contains the shifts of origin that define the multi- */
/* > shift QR sweep. On output S may be reordered. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* > H is COMPLEX array, dimension (LDH,N) */
/* > On input H contains a Hessenberg matrix. On output a */
/* > multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied */
/* > to the isolated diagonal block in rows and columns KTOP */
/* > through KBOT. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > LDH is the leading dimension of H just as declared in the */
/* > calling procedure. LDH >= MAX(1,N). */
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
/* > applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ,IHIZ) */
/* > If WANTZ = .TRUE., then the QR Sweep unitary */
/* > similarity transformation is accumulated into */
/* > Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right. */
/* > If WANTZ = .FALSE., then Z is unreferenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > LDA is the leading dimension of Z just as declared in */
/* > the calling procedure. LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension (LDV,NSHFTS/2) */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > LDV is the leading dimension of V as declared in the */
/* > calling procedure. LDV >= 3. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* > U is COMPLEX array, dimension (LDU,2*NSHFTS) */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > LDU is the leading dimension of U just as declared in the */
/* > in the calling subroutine. LDU >= 2*NSHFTS. */
/* > \endverbatim */
/* > */
/* > \param[in] NV */
/* > \verbatim */
/* > NV is INTEGER */
/* > NV is the number of rows in WV agailable for workspace. */
/* > NV >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WV */
/* > \verbatim */
/* > WV is COMPLEX array, dimension (LDWV,2*NSHFTS) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWV */
/* > \verbatim */
/* > LDWV is INTEGER */
/* > LDWV is the leading dimension of WV as declared in the */
/* > in the calling subroutine. LDWV >= NV. */
/* > \endverbatim */
/* > \param[in] NH */
/* > \verbatim */
/* > NH is INTEGER */
/* > NH is the number of columns in array WH available for */
/* > workspace. NH >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WH */
/* > \verbatim */
/* > WH is COMPLEX array, dimension (LDWH,NH) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWH */
/* > \verbatim */
/* > LDWH is INTEGER */
/* > Leading dimension of WH just as declared in the */
/* > calling procedure. LDWH >= 2*NSHFTS. */
/* > \endverbatim */
/* > */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complexOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Karen Braman and Ralph Byers, Department of Mathematics, */
/* > University of Kansas, USA */
/* > */
/* > Lars Karlsson, Daniel Kressner, and Bruno Lang */
/* > */
/* > Thijs Steel, Department of Computer science, */
/* > KU Leuven, Belgium */
/* > \par References: */
/* ================ */
/* > */
/* > K. Braman, R. Byers and R. Mathias, The Multi-Shift QR */
/* > Algorithm Part I: Maintaining Well Focused Shifts, and Level 3 */
/* > Performance, SIAM Journal of Matrix Analysis, volume 23, pages */
/* > 929--947, 2002. */
/* > */
/* > Lars Karlsson, Daniel Kressner, and Bruno Lang, Optimally packed */
/* > chains of bulges in multishift QR algorithms. */
/* > ACM Trans. Math. Softw. 40, 2, Article 12 (February 2014). */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void claqr5_(logical *wantt, logical *wantz, aocl_int_t *kacc22, aocl_int_t *n, aocl_int_t *ktop,
             aocl_int_t *kbot, aocl_int_t *nshfts, scomplex *s, scomplex *h__, aocl_int_t *ldh,
             aocl_int_t *iloz, aocl_int_t *ihiz, scomplex *z__, aocl_int_t *ldz, scomplex *v,
             aocl_int_t *ldv, scomplex *u, aocl_int_t *ldu, aocl_int_t *nv, scomplex *wv,
             aocl_int_t *ldwv, aocl_int_t *nh, scomplex *wh, aocl_int_t *ldwh)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_claqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h__, ldh, iloz, ihiz, z__,
                       ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
#else
    aocl_int64_t kacc22_64 = *kacc22;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ktop_64 = *ktop;
    aocl_int64_t kbot_64 = *kbot;
    aocl_int64_t nshfts_64 = *nshfts;
    aocl_int64_t ldh_64 = *ldh;
    aocl_int64_t iloz_64 = *iloz;
    aocl_int64_t ihiz_64 = *ihiz;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t ldv_64 = *ldv;
    aocl_int64_t ldu_64 = *ldu;
    aocl_int64_t nv_64 = *nv;
    aocl_int64_t ldwv_64 = *ldwv;
    aocl_int64_t nh_64 = *nh;
    aocl_int64_t ldwh_64 = *ldwh;

    aocl_lapack_claqr5(wantt, wantz, &kacc22_64, &n_64, &ktop_64, &kbot_64, &nshfts_64, s, h__,
                       &ldh_64, &iloz_64, &ihiz_64, z__, &ldz_64, v, &ldv_64, u, &ldu_64, &nv_64,
                       wv, &ldwv_64, &nh_64, wh, &ldwh_64);
#endif
}

void aocl_lapack_claqr5(logical *wantt, logical *wantz, aocl_int64_t *kacc22, aocl_int64_t *n,
                        aocl_int64_t *ktop, aocl_int64_t *kbot, aocl_int64_t *nshfts, scomplex *s,
                        scomplex *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                        scomplex *z__, aocl_int64_t *ldz, scomplex *v, aocl_int64_t *ldv, scomplex *u,
                        aocl_int64_t *ldu, aocl_int64_t *nv, scomplex *wv, aocl_int64_t *ldwv,
                        aocl_int64_t *nh, scomplex *wh, aocl_int64_t *ldwh)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(
        buffer, 256,
        "claqr5 inputs: kacc22 %lld, n %lld, ktop %lld, kbot %lld, nshfts %lld, ldh %lld, iloz "
        "%lld, ihiz %lld, ldz %lld, ldv %lld, ldu %lld, nv %lld, ldwv %lld, nh %lld, ldwh %lld",
        *kacc22, *n, *ktop, *kbot, *nshfts, *ldh, *iloz, *ihiz, *ldz, *ldv, *ldu, *nv, *ldwv, *nh,
        *ldwh);
#else
    snprintf(buffer, 256,
             "claqr5 inputs: kacc22 %d, n %d, ktop %d, kbot %d, nshfts %d, ldh %d, iloz %d, ihiz "
             "%d, ldz %d, ldv %d, ldu %d, nv %d, ldwv %d, nh %d, ldwh %d",
             *kacc22, *n, *ktop, *kbot, *nshfts, *ldh, *iloz, *ihiz, *ldz, *ldv, *ldu, *nv, *ldwv,
             *nh, *ldwh);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t h_dim1, h_offset, u_dim1, u_offset, v_dim1, v_offset, wh_dim1, wh_offset, wv_dim1,
        wv_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10,
        i__11;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10;
    scomplex q__1, q__2, q__3, q__4, q__5, q__6, q__7;
    /* Local variables */
    aocl_int64_t j, k, m, i2, k1, i4;
    real h11, h12, h21, h22;
    aocl_int64_t m22, ns, nu;
    scomplex vt[3];
    real scl;
    aocl_int64_t kdu, kms;
    real ulp, tst1, tst2;
    scomplex beta;
    logical bmp22;
    aocl_int64_t jcol, jlen, jbot, mbot, jtop, jrow, mtop;
    scomplex alpha;
    logical accum;
    aocl_int64_t ndcol, incol, krcol, nbmps;
    extern real slamch_(char *);
    real safmin;
    scomplex refsum;
    real smlnum;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* ==== If there are no shifts, then there is nothing to do. ==== */
    /* Parameter adjustments */
    --s;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    wv_dim1 = *ldwv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    wh_dim1 = *ldwh;
    wh_offset = 1 + wh_dim1;
    wh -= wh_offset;
    /* Function Body */
    q__5.real = 0.f;
    q__5.imag = 0.f;
    if(*nshfts < 2)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* ==== If the active block is empty or 1-by-1, then there */
    /* . is nothing to do. ==== */
    if(*ktop >= *kbot)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* ==== NSHFTS is supposed to be even, but if it is odd, */
    /* . then simply reduce it by one. ==== */
    ns = *nshfts - *nshfts % 2;
    /* ==== Machine constants for deflation ==== */
    safmin = slamch_("SAFE MINIMUM");
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real)(*n) / ulp);
    /* ==== Use accumulated reflections to update far-from-diagonal */
    /* . entries ? ==== */
    accum = *kacc22 == 1 || *kacc22 == 2;
    /* ==== clear trash ==== */
    if(*ktop + 2 <= *kbot)
    {
        i__1 = *ktop + 2 + *ktop * h_dim1;
        h__[i__1].real = 0.f;
        h__[i__1].imag = 0.f; // , expr subst
    }
    /* ==== NBMPS = number of 2-shift bulges in the chain ==== */
    nbmps = ns / 2;
    /* ==== KDU = width of slab ==== */
    kdu = nbmps << 2;
    /* ==== Create and chase chains of NBMPS bulges ==== */
    i__1 = *kbot - 2;
    i__2 = nbmps << 1;
    for(incol = *ktop - (nbmps << 1) + 1; i__2 < 0 ? incol >= i__1 : incol <= i__1; incol += i__2)
    {
        /* JTOP = Index from which updates from the right start. */
        if(accum)
        {
            jtop = fla_max(*ktop, incol);
        }
        else if(*wantt)
        {
            jtop = 1;
        }
        else
        {
            jtop = *ktop;
        }
        ndcol = incol + kdu;
        if(accum)
        {
            aocl_lapack_claset("ALL", &kdu, &kdu, &c_b1, &c_b2, &u[u_offset], ldu);
        }
        /* ==== Near-the-diagonal bulge chase. The following loop */
        /* . performs the near-the-diagonal part of a small bulge */
        /* . multi-shift QR sweep. Each 4*NBMPS column diagonal */
        /* . chunk extends from column INCOL to column NDCOL */
        /* . (including both column INCOL and column NDCOL). The */
        /* . following loop chases a 2*NBMPS+1 column long chain of */
        /* . NBMPS bulges 2*NBMPS columns to the right. (INCOL */
        /* . may be less than KTOP and and NDCOL may be greater than */
        /* . KBOT indicating phantom columns from which to chase */
        /* . bulges before they are actually introduced or to which */
        /* . to chase bulges beyond column KBOT.) ==== */
        /* Computing MIN */
        i__4 = incol + (nbmps << 1) - 1;
        i__5 = *kbot - 2; // , expr subst
        i__3 = fla_min(i__4, i__5);
        for(krcol = incol; krcol <= i__3; ++krcol)
        {
            /* ==== Bulges number MTOP to MBOT are active double implicit */
            /* . shift bulges. There may or may not also be small */
            /* . 2-by-2 bulge, if there is room. The inactive bulges */
            /* . (if any) must wait until the active bulges have moved */
            /* . down the diagonal to make room. The phantom matrix */
            /* . paradigm described above helps keep track. ==== */
            /* Computing MAX */
            real v1r, v1i, v2r, v2i, v3r, v3i;
            real u1r, u1i, u2r, u2i, u3r, u3i, u4i, u4r;

            i__4 = 1;
            i__5 = (*ktop - krcol) / 2 + 1; // , expr subst
            mtop = fla_max(i__4, i__5);
            /* Computing MIN */
            i__4 = nbmps;
            i__5 = (*kbot - krcol - 1) / 2; // , expr subst
            mbot = fla_min(i__4, i__5);
            m22 = mbot + 1;
            bmp22 = mbot < nbmps && krcol + (m22 - 1 << 1) == *kbot - 2;
            /* ==== Generate reflections to chase the chain right */
            /* . one column. (The minimum value of K is KTOP-1.) ==== */
            if(bmp22)
            {
                /* ==== Special case: 2-by-2 reflection at bottom treated */
                /* . separately ==== */
                k = krcol + (m22 - 1 << 1);
                if(k == *ktop - 1)
                {
                    i__4 = m22 * v_dim1 + 1;
                    aocl_lapack_claqr1(&c__2, &h__[k + 1 + (k + 1) * h_dim1], ldh,
                                       &s[(m22 << 1) - 1], &s[m22 * 2], &v[i__4]);
                    beta.real = v[i__4].real;
                    beta.imag = v[i__4].imag; // , expr subst
                    aocl_lapack_clarfg(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[i__4]);
                }
                else
                {
                    i__4 = k + 1 + k * h_dim1;
                    beta.real = h__[i__4].real;
                    beta.imag = h__[i__4].imag; // , expr subst
                    i__4 = m22 * v_dim1 + 2;
                    i__5 = k + 2 + k * h_dim1;
                    v[i__4].real = h__[i__5].real;
                    v[i__4].imag = h__[i__5].imag; // , expr subst
                    aocl_lapack_clarfg(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1,
                                       &v[m22 * v_dim1 + 1]);
                    i__4 = k + 1 + k * h_dim1;
                    h__[i__4].real = beta.real;
                    h__[i__4].imag = beta.imag; // , expr subst
                    i__4 = k + 2 + k * h_dim1;
                    h__[i__4].real = 0.f;
                    h__[i__4].imag = 0.f; // , expr subst
                }
                /* ==== Perform update from right within */
                /* . computational window. ==== */
                /* Computing MIN */
                i__5 = *kbot;
                i__6 = k + 3; // , expr subst
                i__4 = fla_min(i__5, i__6);
                v1r = v[m22 * v_dim1 + 1].real;
                v1i = v[m22 * v_dim1 + 1].imag;
                v2r = v[m22 * v_dim1 + 2].real;
                v2i = v[m22 * v_dim1 + 2].imag;
                v3r = v[m22 * v_dim1 + 3].real;
                v3i = v[m22 * v_dim1 + 3].imag;

                for(j = jtop; j <= i__4; ++j)
                {
                    i__6 = j + (k + 1) * h_dim1;
                    i__8 = j + (k + 2) * h_dim1;
                    u1r = h__[i__6].real;
                    u1i = h__[i__6].imag;
                    u2r = h__[i__8].real;
                    u2i = h__[i__8].imag;
                    q__3.real = v2r * u2r - v2i * u2i;
                    q__3.imag = v2r * u2i + v2i * u2r; // , expr subst
                    q__2.real = u1r + q__3.real;
                    q__2.imag = u1i + q__3.imag; // , expr subst
                    q__1.real = v1r * q__2.real - v1i * q__2.imag;
                    q__1.imag = v1r * q__2.imag + v1i * q__2.real; // , expr subst
                    refsum.real = q__1.real;
                    refsum.imag = q__1.imag; // , expr subst
                    q__1.real = u1r - refsum.real;
                    q__1.imag = u1i - refsum.imag; // , expr subst
                    h__[i__6].real = q__1.real;
                    h__[i__6].imag = q__1.imag; // , expr subst
                    q__2.real = refsum.real * v2r + refsum.imag * v2i;
                    q__2.imag = -refsum.real * v2i + refsum.imag * v2r; // , expr subst
                    q__1.real = u2r - q__2.real;
                    q__1.imag = u2i - q__2.imag; // , expr subst
                    h__[i__8].real = q__1.real;
                    h__[i__8].imag = q__1.imag; // , expr subst
                    /* L30: */
                }
                /* ==== Perform update from left within */
                /* . computational window. ==== */
                if(accum)
                {
                    jbot = fla_min(ndcol, *kbot);
                }
                else if(*wantt)
                {
                    jbot = *n;
                }
                else
                {
                    jbot = *kbot;
                }
                i__4 = jbot;
                for(j = k + 1; j <= i__4; ++j)
                {
                    i__5 = k + 1 + j * h_dim1;
                    i__6 = k + 2 + j * h_dim1;
                    u1r = h__[i__5].real;
                    u1i = h__[i__5].imag;
                    u2r = h__[i__6].real;
                    u2i = h__[i__6].imag;
                    q__4.real = v2r * u2r + v2i * u2i;
                    q__4.imag = v2r * u2i - v2i * u2r; // , expr subst
                    q__3.real = u1r + q__4.real;
                    q__3.imag = u1i + q__4.imag; // , expr subst
                    q__1.real = v1r * q__3.real + v1i * q__3.imag;
                    q__1.imag = v1r * q__3.imag - v1i * q__3.real; // , expr subst
                    refsum.real = q__1.real;
                    refsum.imag = q__1.imag; // , expr subst
                    q__1.real = u1r - refsum.real;
                    q__1.imag = u1i - refsum.imag; // , expr subst
                    h__[i__5].real = q__1.real;
                    h__[i__5].imag = q__1.imag; // , expr subst
                    q__2.real = refsum.real * v2r - refsum.imag * v2i;
                    q__2.imag = refsum.real * v2i + refsum.imag * v2r; // , expr subst
                    q__1.real = u2r - q__2.real;
                    q__1.imag = u2i - q__2.imag; // , expr subst
                    h__[i__6].real = q__1.real;
                    h__[i__6].imag = q__1.imag; // , expr subst
                    /* L40: */
                }
                /* ==== The following convergence test requires that */
                /* . the tradition small-compared-to-nearby-diagonals */
                /* . criterion and the Ahues & Tisseur (LAWN 122, 1997) */
                /* . criteria both be satisfied. The latter improves */
                /* . accuracy in some examples. Falling back on an */
                /* . alternate convergence criterion when TST1 or TST2 */
                /* . is zero (as done here) is traditional but probably */
                /* . unnecessary. ==== */
                if(k >= *ktop)
                {
                    i__4 = k + 1 + k * h_dim1;
                    u1r = h__[i__4].real;
                    u1i = h__[i__4].imag;
                    u2r = h__[k + (k + 1) * h_dim1].real;
                    u2i = h__[k + (k + 1) * h_dim1].imag;
                    u3r = h__[k + k * h_dim1].real;
                    u3i = h__[k + k * h_dim1].imag;
                    u4r = h__[k + 1 + (k + 1) * h_dim1].real;
                    u4i = h__[k + 1 + (k + 1) * h_dim1].imag;
                    if(u1r != 0.f || u1i != 0.f)
                    {
                        tst1 = (r__1 = u3r, f2c_abs(r__1)) + (r__2 = u3i, f2c_abs(r__2))
                               + ((r__3 = u4r, f2c_abs(r__3)) + (r__4 = u4i, f2c_abs(r__4)));
                        if(tst1 == 0.f)
                        {
                            if(k >= *ktop + 1)
                            {
                                i__4 = k + (k - 1) * h_dim1;
                                tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                        + (r__2 = h__[k + (k - 1) * h_dim1].imag, f2c_abs(r__2));
                            }
                            if(k >= *ktop + 2)
                            {
                                i__4 = k + (k - 2) * h_dim1;
                                tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                        + (r__2 = h__[k + (k - 2) * h_dim1].imag, f2c_abs(r__2));
                            }
                            if(k >= *ktop + 3)
                            {
                                i__4 = k + (k - 3) * h_dim1;
                                tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                        + (r__2 = h__[k + (k - 3) * h_dim1].imag, f2c_abs(r__2));
                            }
                            if(k <= *kbot - 2)
                            {
                                i__4 = k + 2 + (k + 1) * h_dim1;
                                tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                        + (r__2 = h__[k + 2 + (k + 1) * h_dim1].imag, f2c_abs(r__2));
                            }
                            if(k <= *kbot - 3)
                            {
                                i__4 = k + 3 + (k + 1) * h_dim1;
                                tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                        + (r__2 = h__[k + 3 + (k + 1) * h_dim1].imag, f2c_abs(r__2));
                            }
                            if(k <= *kbot - 4)
                            {
                                i__4 = k + 4 + (k + 1) * h_dim1;
                                tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                        + (r__2 = h__[k + 4 + (k + 1) * h_dim1].imag, f2c_abs(r__2));
                            }
                        }
                        /* Computing MAX */
                        r__3 = smlnum;
                        r__4 = ulp * tst1; // , expr subst
                        if((r__1 = u1r, f2c_abs(r__1)) + (r__2 = u1i, f2c_abs(r__2))
                           <= fla_max(r__3, r__4))
                        {
                            /* Computing MAX */
                            r__5 = (r__1 = u1r, f2c_abs(r__1)) + (r__2 = u1i, f2c_abs(r__2));
                            r__6 = (r__3 = u2r, f2c_abs(r__3))
                                   + (r__4 = u2i, f2c_abs(r__4)); // , expr subst
                            h12 = fla_max(r__5, r__6);
                            /* Computing MIN */
                            r__5 = (r__1 = u1r, f2c_abs(r__1)) + (r__2 = u1i, f2c_abs(r__2));
                            r__6 = (r__3 = u2r, f2c_abs(r__3))
                                   + (r__4 = u2i, f2c_abs(r__4)); // , expr subst
                            h21 = fla_min(r__5, r__6);
                            q__2.real = u3r - u4r;
                            q__2.imag = u3i - u4i; // , expr subst
                            q__1.real = q__2.real;
                            q__1.imag = q__2.imag; // , expr subst
                            /* Computing MAX */
                            r__5 = (r__1 = u4r, f2c_abs(r__1)) + (r__2 = u4i, f2c_abs(r__2));
                            r__6 = (r__3 = q__1.real, f2c_abs(r__3))
                                   + (r__4 = q__1.imag, f2c_abs(r__4)); // , expr subst
                            h11 = fla_max(r__5, r__6);
                            q__2.real = u3r - u4r;
                            q__2.imag = u3i - u4i; // , expr subst
                            q__1.real = q__2.real;
                            q__1.imag = q__2.imag; // , expr subst
                            /* Computing MIN */
                            r__5 = (r__1 = u4r, f2c_abs(r__1)) + (r__2 = u4i, f2c_abs(r__2));
                            r__6 = (r__3 = q__1.real, f2c_abs(r__3))
                                   + (r__4 = q__1.imag, f2c_abs(r__4)); // , expr subst
                            h22 = fla_min(r__5, r__6);
                            scl = h11 + h12;
                            tst2 = h22 * (h11 / scl);
                            /* Computing MAX */
                            r__1 = smlnum;
                            r__2 = ulp * tst2; // , expr subst
                            if(tst2 == 0.f || h21 * (h12 / scl) <= fla_max(r__1, r__2))
                            {
                                i__4 = k + 1 + k * h_dim1;
                                h__[i__4].real = 0.f;
                                h__[i__4].imag = 0.f; // , expr subst
                            }
                        }
                    }
                }
                /* ==== Accumulate orthogonal transformations. ==== */
                if(accum)
                {
                    kms = k - incol;
                    /* Computing MAX */
                    i__4 = 1;
                    i__5 = *ktop - incol; // , expr subst
                    i__6 = kdu;

                    v1r = v[m22 * v_dim1 + 1].real;
                    v1i = v[m22 * v_dim1 + 1].imag;
                    v2r = v[m22 * v_dim1 + 2].real;
                    v2i = v[m22 * v_dim1 + 2].imag;

                    for(j = fla_max(i__4, i__5); j <= i__6; ++j)
                    {
                        i__5 = j + (kms + 1) * u_dim1;
                        i__8 = j + (kms + 2) * u_dim1;

                        u1r = u[i__5].real;
                        u1i = u[i__5].imag;
                        u2r = u[i__8].real;
                        u2i = u[i__8].imag;
                        q__3.real = v2r * u2r - v2i * u2i;
                        q__3.imag = v2r * u2i + v2i * u2r; // , expr subst
                        q__2.real = u1r + q__3.real;
                        q__2.imag = u1i + q__3.imag; // , expr subst
                        q__1.real = v1r * q__2.real - v1i * q__2.imag;
                        q__1.imag = v1r * q__2.imag + v1i * q__2.real; // , expr subst
                        refsum.real = q__1.real;
                        refsum.imag = q__1.imag; // , expr subst
                        q__1.real = u1r - refsum.real;
                        q__1.imag = u1i - refsum.imag; // , expr subst
                        u[i__5].real = q__1.real;
                        u[i__5].imag = q__1.imag; // , expr subst
                        q__2.real = refsum.real * v2r + refsum.imag * v2i;
                        q__2.imag = -refsum.real * v2i + refsum.imag * v2r; // , expr subst
                        q__1.real = u2r - q__2.real;
                        q__1.imag = u2i - q__2.imag; // , expr subst
                        u[i__8].real = q__1.real;
                        u[i__8].imag = q__1.imag; // , expr subst
                        /* L50: */
                    }
                }
                else if(*wantz)
                {
                    i__6 = *ihiz;
                    for(j = *iloz; j <= i__6; ++j)
                    {
                        i__4 = m22 * v_dim1 + 1;
                        i__5 = j + (k + 1) * z_dim1;
                        i__7 = m22 * v_dim1 + 2;
                        i__8 = j + (k + 2) * z_dim1;
                        q__3.real = v[i__7].real * z__[i__8].real - v[i__7].imag * z__[i__8].imag;
                        q__3.imag = v[i__7].real * z__[i__8].imag + v[i__7].imag * z__[i__8].real; // , expr subst
                        q__2.real = z__[i__5].real + q__3.real;
                        q__2.imag = z__[i__5].imag + q__3.imag; // , expr subst
                        q__1.real = v[i__4].real * q__2.real - v[i__4].imag * q__2.imag;
                        q__1.imag = v[i__4].real * q__2.imag + v[i__4].imag * q__2.real; // , expr subst
                        refsum.real = q__1.real;
                        refsum.imag = q__1.imag; // , expr subst
                        i__4 = j + (k + 1) * z_dim1;
                        i__5 = j + (k + 1) * z_dim1;
                        q__1.real = z__[i__5].real - refsum.real;
                        q__1.imag = z__[i__5].imag - refsum.imag; // , expr subst
                        z__[i__4].real = q__1.real;
                        z__[i__4].imag = q__1.imag; // , expr subst
                        i__4 = j + (k + 2) * z_dim1;
                        i__5 = j + (k + 2) * z_dim1;
                        q__3.real = v[m22 * v_dim1 + 2].real;
                        q__3.imag = -v[m22 * v_dim1 + 2].imag;
                        q__2.real = refsum.real * q__3.real - refsum.imag * q__3.imag;
                        q__2.imag = refsum.real * q__3.imag + refsum.imag * q__3.real; // , expr subst
                        q__1.real = z__[i__5].real - q__2.real;
                        q__1.imag = z__[i__5].imag - q__2.imag; // , expr subst
                        z__[i__4].real = q__1.real;
                        z__[i__4].imag = q__1.imag; // , expr subst
                        /* L60: */
                    }
                }
            }
            /* ==== Normal case: Chain of 3-by-3 reflections ==== */
            i__6 = mtop;
            for(m = mbot; m >= i__6; --m)
            {

                v1r = v[m * v_dim1 + 1].real;
                v1i = v[m * v_dim1 + 1].imag;
                v2r = v[m * v_dim1 + 2].real;
                v2i = v[m * v_dim1 + 2].imag;
                v3r = v[m * v_dim1 + 3].real;
                v3i = v[m * v_dim1 + 3].imag;

                k = krcol + (m - 1 << 1);
                if(k == *ktop - 1)
                {
                    aocl_lapack_claqr1(&c__3, &h__[*ktop + *ktop * h_dim1], ldh, &s[(m << 1) - 1],
                                       &s[m * 2], &v[m * v_dim1 + 1]);
                    i__4 = m * v_dim1 + 1;
                    alpha.real = v[i__4].real;
                    alpha.imag = v[i__4].imag; // , expr subst
                    aocl_lapack_clarfg(&c__3, &alpha, &v[m * v_dim1 + 2], &c__1,
                                       &v[m * v_dim1 + 1]);
                }
                else
                {
                    /* ==== Perform delayed transformation of row below */
                    /* . Mth bulge. Exploit fact that first two elements */
                    /* . of row are actually zero. ==== */
                    q__2.real = v1r * v3r - v1i * v3i;
                    q__2.imag = v1r * v3i + v1i * v3r; // , expr subst
                    i__4 = k + 3 + (k + 1) * h_dim1;
                    i__7 = k + 3 + (k + 2) * h_dim1;
                    u1r = h__[i__4].real;
                    u1i = h__[i__4].imag;
                    u2r = h__[i__7].real;
                    u2i = h__[i__7].imag;
                    q__1.real = q__2.real * u2r - q__2.imag * u2i;
                    q__1.imag = q__2.real * u2i + q__2.imag * u2r; // , expr subst
                    refsum.real = q__1.real;
                    refsum.imag = q__1.imag; // , expr subst
                    i__4 = k + 3 + k * h_dim1;
                    q__1.real = -refsum.real;
                    q__1.imag = -refsum.imag; // , expr subst
                    h__[i__4].real = q__1.real;
                    h__[i__4].imag = q__1.imag; // , expr subst
                    i__4 = k + 3 + (k + 1) * h_dim1;
                    q__2.real = -refsum.real;
                    q__2.imag = -refsum.imag; // , expr subst
                    q__1.real = q__2.real * v2r + q__2.imag * v2i;
                    q__1.imag = -q__2.real * v2i + q__2.imag * v2r; // , expr subst
                    h__[i__4].real = q__1.real;
                    h__[i__4].imag = q__1.imag; // , expr subst
                    u1r = q__1.real;
                    u1i = q__1.imag;
                    q__2.real = refsum.real * v3r + refsum.imag * v3i;
                    q__2.imag = -refsum.real * v3i + refsum.imag * v3r; // , expr subst
                    q__1.real = u2r - q__2.real;
                    q__1.imag = u2i - q__2.imag; // , expr subst
                    h__[i__7].real = q__1.real;
                    h__[i__7].imag = q__1.imag; // , expr subst
                    u2r = q__1.real;
                    u2i = q__1.imag;
                    /* ==== Calculate reflection to move */
                    /* . Mth bulge one step. ==== */
                    i__4 = k + 1 + k * h_dim1;
                    beta.real = h__[i__4].real;
                    beta.imag = h__[i__4].imag; // , expr subst
                    i__4 = m * v_dim1 + 2;
                    i__5 = k + 2 + k * h_dim1;
                    v[i__4].real = h__[i__5].real;
                    v[i__4].imag = h__[i__5].imag; // , expr subst
                    i__4 = m * v_dim1 + 3;
                    i__5 = k + 3 + k * h_dim1;
                    u3r = h__[i__5].real;
                    u3i = h__[i__5].imag;
                    v[i__4].real = u3r;
                    v[i__4].imag = u3i; // , expr subst
                    aocl_lapack_clarfg(&c__3, &beta, &v[m * v_dim1 + 2], &c__1, &v[m * v_dim1 + 1]);
                    /* ==== A Bulge may collapse because of vigilant */
                    /* . deflation or destructive underflow. In the */
                    /* . underflow case, try the two-small-subdiagonals */
                    /* . trick to try to reinflate the bulge. ==== */
                    if(u3r != 0.f || u3i != 0.f || (u1r != 0.f || u1i != 0.f)
                       || u2r == 0.f && u2i == 0.f)
                    {
                        /* ==== Typical case: not collapsed (yet). ==== */
                        i__4 = k + 1 + k * h_dim1;
                        h__[i__4].real = beta.real;
                        h__[i__4].imag = beta.imag; // , expr subst
                        i__4 = k + 2 + k * h_dim1;
                        h__[i__4].real = 0.f;
                        h__[i__4].imag = 0.f; // , expr subst
                        i__4 = k + 3 + k * h_dim1;
                        h__[i__4].real = 0.f;
                        h__[i__4].imag = 0.f; // , expr subst
                    }
                    else
                    {
                        /* ==== Atypical case: collapsed. Attempt to */
                        /* . reintroduce ignoring H(K+1,K) and H(K+2,K). */
                        /* . If the fill resulting from the new */
                        /* . reflector is too large, then abandon it. */
                        /* . Otherwise, use the new one. ==== */
                        aocl_lapack_claqr1(&c__3, &h__[k + 1 + (k + 1) * h_dim1], ldh,
                                           &s[(m << 1) - 1], &s[m * 2], vt);
                        alpha.real = vt[0].real;
                        alpha.imag = vt[0].imag; // , expr subst
                        aocl_lapack_clarfg(&c__3, &alpha, &vt[1], &c__1, vt);
                        q__2.real = vt->real;
                        q__2.imag = -vt->imag;
                        i__4 = k + 1 + k * h_dim1;
                        q__3.real = vt[1].real;
                        q__3.imag = -vt[1].imag;
                        i__5 = k + 2 + k * h_dim1;
                        q__4.real = q__5.real * h__[i__5].real - q__5.imag * h__[i__5].imag;
                        q__4.imag = q__5.real * h__[i__5].imag + q__5.imag * h__[i__5].real; // , expr subst
                        q__3.real = h__[i__4].real + q__4.real;
                        q__3.imag = h__[i__4].imag + q__4.imag; // , expr subst
                        q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                        q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                        refsum.real = q__1.real;
                        refsum.imag = q__1.imag; // , expr subst
                        i__4 = k + 2 + k * h_dim1;
                        q__3.real = refsum.real * vt[1].real - refsum.imag * vt[1].imag;
                        q__3.imag = refsum.real * vt[1].imag + refsum.imag * vt[1].real; // , expr subst
                        q__2.real = h__[i__4].real - q__3.real;
                        q__2.imag = h__[i__4].imag - q__3.imag; // , expr subst
                        q__1.real = q__2.real;
                        q__1.imag = q__2.imag; // , expr subst
                        q__5.real = refsum.real * vt[2].real - refsum.imag * vt[2].imag;
                        q__5.imag = refsum.real * vt[2].imag + refsum.imag * vt[2].real; // , expr subst
                        q__4.real = q__5.real;
                        q__4.imag = q__5.imag; // , expr subst
                        i__5 = k + k * h_dim1;
                        i__7 = k + 1 + (k + 1) * h_dim1;
                        i__8 = k + 2 + (k + 2) * h_dim1;
                        if((r__1 = q__1.real, f2c_abs(r__1)) + (r__2 = q__1.imag, f2c_abs(r__2))
                               + ((r__3 = q__4.real, f2c_abs(r__3)) + (r__4 = q__4.imag, f2c_abs(r__4)))
                           > ulp
                                 * ((r__5 = h__[i__5].real, f2c_abs(r__5))
                                    + (r__6 = h__[k + k * h_dim1].imag, f2c_abs(r__6))
                                    + ((r__7 = h__[i__7].real, f2c_abs(r__7))
                                       + (r__8 = h__[k + 1 + (k + 1) * h_dim1].imag, f2c_abs(r__8)))
                                    + ((r__9 = h__[i__8].real, f2c_abs(r__9))
                                       + (r__10 = h__[k + 2 + (k + 2) * h_dim1].imag,
                                          f2c_abs(r__10)))))
                        {
                            /* ==== Starting a new bulge here would */
                            /* . create non-negligible fill. Use */
                            /* . the old one with trepidation. ==== */
                            i__4 = k + 1 + k * h_dim1;
                            h__[i__4].real = beta.real;
                            h__[i__4].imag = beta.imag; // , expr subst
                            i__4 = k + 2 + k * h_dim1;
                            h__[i__4].real = 0.f;
                            h__[i__4].imag = 0.f; // , expr subst
                            i__4 = k + 3 + k * h_dim1;
                            h__[i__4].real = 0.f;
                            h__[i__4].imag = 0.f; // , expr subst
                        }
                        else
                        {
                            /* ==== Starting a new bulge here would */
                            /* . create only negligible fill. */
                            /* . Replace the old reflector with */
                            /* . the new one. ==== */
                            i__4 = k + 1 + k * h_dim1;
                            i__5 = k + 1 + k * h_dim1;
                            q__1.real = h__[i__5].real - refsum.real;
                            q__1.imag = h__[i__5].imag - refsum.imag; // , expr subst
                            h__[i__4].real = q__1.real;
                            h__[i__4].imag = q__1.imag; // , expr subst
                            i__4 = k + 2 + k * h_dim1;
                            h__[i__4].real = 0.f;
                            h__[i__4].imag = 0.f; // , expr subst
                            i__4 = k + 3 + k * h_dim1;
                            h__[i__4].real = 0.f;
                            h__[i__4].imag = 0.f; // , expr subst
                            i__4 = m * v_dim1 + 1;
                            v[i__4].real = vt[0].real;
                            v[i__4].imag = vt[0].imag; // , expr subst
                            i__4 = m * v_dim1 + 2;
                            v[i__4].real = vt[1].real;
                            v[i__4].imag = vt[1].imag; // , expr subst
                            i__4 = m * v_dim1 + 3;
                            v[i__4].real = vt[2].real;
                            v[i__4].imag = vt[2].imag; // , expr subst
                        }
                    }
                }
                /* ==== Apply reflection from the right and */
                /* . the first column of update from the left. */
                /* . These updates are required for the vigilant */
                /* . deflation check. We still delay most of the */
                /* . updates from the left for efficiency. ==== */
                /* Computing MIN */
                i__5 = *kbot;
                i__7 = k + 3; // , expr subst
                i__4 = fla_min(i__5, i__7);
                v1r = v[m * v_dim1 + 1].real;
                v1i = v[m * v_dim1 + 1].imag;
                v2r = v[m * v_dim1 + 2].real;
                v2i = v[m * v_dim1 + 2].imag;
                v3r = v[m * v_dim1 + 3].real;
                v3i = v[m * v_dim1 + 3].imag;
                for(j = jtop; j <= i__4; ++j)
                {
                    i__7 = j + (k + 1) * h_dim1;
                    i__9 = j + (k + 2) * h_dim1;
                    i__11 = j + (k + 3) * h_dim1;

                    u1r = h__[i__7].real;
                    u1i = h__[i__7].imag;
                    u2r = h__[i__9].real;
                    u2i = h__[i__9].imag;
                    u3r = h__[i__11].real;
                    u3i = h__[i__11].imag;

                    q__4.real = v2r * u2r - v2i * u2i;
                    q__4.imag = v2r * u2i + v2i * u2r; // , expr subst
                    q__3.real = u1r + q__4.real;
                    q__3.imag = u1i + q__4.imag; // , expr subst

                    q__5.real = v3r * u3r - v3i * u3i;
                    q__5.imag = v3r * u3i + v3i * u3r; // , expr subst
                    q__2.real = q__3.real + q__5.real;
                    q__2.imag = q__3.imag + q__5.imag; // , expr subst
                    q__1.real = v1r * q__2.real - v1i * q__2.imag;
                    q__1.imag = v1r * q__2.imag + v1i * q__2.real; // , expr subst
                    refsum.real = q__1.real;
                    refsum.imag = q__1.imag; // , expr subst
                    q__1.real = u1r - refsum.real;
                    q__1.imag = u1i - refsum.imag; // , expr subst
                    h__[i__7].real = q__1.real;
                    h__[i__7].imag = q__1.imag; // , expr subst
                    q__2.real = refsum.real * v2r + refsum.imag * v2i;
                    q__2.imag = -refsum.real * v2i + refsum.imag * v2r; // , expr subst
                    q__1.real = u2r - q__2.real;
                    q__1.imag = u2i - q__2.imag; // , expr subst
                    h__[i__9].real = q__1.real;
                    h__[i__9].imag = q__1.imag; // , expr subst
                    q__2.real = refsum.real * v3r + refsum.imag * v3i;
                    q__2.imag = -refsum.real * v3i + refsum.imag * v3r; // , expr subst
                    q__1.real = u3r - q__2.real;
                    q__1.imag = u3i - q__2.imag; // , expr subst
                    h__[i__11].real = q__1.real;
                    h__[i__11].imag = q__1.imag; // , expr subst
                    /* L70: */
                }
                /* ==== Perform update from left for subsequent */
                /* . column. ==== */
                q__2.real = v1r;
                q__2.imag = -v1i;
                q__6.real = v2r;
                q__6.imag = -v2i;
                i__4 = k + 1 + (k + 1) * h_dim1;
                i__5 = k + 2 + (k + 1) * h_dim1;
                i__11 = k + 3 + (k + 1) * h_dim1;
                u1r = h__[i__4].real;
                u1i = h__[i__4].imag;
                u2r = h__[i__5].real;
                u2i = h__[i__5].imag;
                u3r = h__[i__11].real;
                u3i = h__[i__11].imag;
                q__5.real = q__6.real * u2r - q__6.imag * u2i;
                q__5.imag = q__6.real * u2i + q__6.imag * u2r; // , expr subst
                q__4.real = u1r + q__5.real;
                q__4.imag = u1i + q__5.imag; // , expr subst
                q__7.real = v3r * u3r + v3i * u3i;
                q__7.imag = v3r * u3i - v3i * u3r; // , expr subst
                q__3.real = q__4.real + q__7.real;
                q__3.imag = q__4.imag + q__7.imag; // , expr subst
                q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                refsum.real = q__1.real;
                refsum.imag = q__1.imag; // , expr subst
                q__1.real = u1r - refsum.real;
                q__1.imag = u1i - refsum.imag; // , expr subst
                h__[i__4].real = q__1.real;
                h__[i__4].imag = q__1.imag; // , expr subst
                q__2.real = refsum.real * v2r - refsum.imag * v2i;
                q__2.imag = refsum.real * v2i + refsum.imag * v2r; // , expr subst
                q__1.real = u2r - q__2.real;
                q__1.imag = u2i - q__2.imag; // , expr subst
                h__[i__5].real = q__1.real;
                h__[i__5].imag = q__1.imag; // , expr subst
                q__2.real = refsum.real * v3r - refsum.imag * v3i;
                q__2.imag = refsum.real * v3i + refsum.imag * v3r; // , expr subst
                q__1.real = u3r - q__2.real;
                q__1.imag = u3i - q__2.imag; // , expr subst
                h__[i__11].real = q__1.real;
                h__[i__11].imag = q__1.imag; // , expr subst
                /* ==== The following convergence test requires that */
                /* . the tradition small-compared-to-nearby-diagonals */
                /* . criterion and the Ahues & Tisseur (LAWN 122, 1997) */
                /* . criteria both be satisfied. The latter improves */
                /* . accuracy in some examples. Falling back on an */
                /* . alternate convergence criterion when TST1 or TST2 */
                /* . is zero (as done here) is traditional but probably */
                /* . unnecessary. ==== */
                if(k < *ktop)
                {
                    continue;
                }
                i__4 = k + 1 + k * h_dim1;
                if(h__[i__4].real != 0.f || h__[i__4].imag != 0.f)
                {
                    i__4 = k + k * h_dim1;
                    i__5 = k + 1 + (k + 1) * h_dim1;
                    tst1 = (r__1 = h__[i__4].real, f2c_abs(r__1))
                           + (r__2 = h__[k + k * h_dim1].imag, f2c_abs(r__2))
                           + ((r__3 = h__[i__5].real, f2c_abs(r__3))
                              + (r__4 = h__[k + 1 + (k + 1) * h_dim1].imag, f2c_abs(r__4)));
                    if(tst1 == 0.f)
                    {
                        if(k >= *ktop + 1)
                        {
                            i__4 = k + (k - 1) * h_dim1;
                            tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                    + (r__2 = h__[k + (k - 1) * h_dim1].imag, f2c_abs(r__2));
                        }
                        if(k >= *ktop + 2)
                        {
                            i__4 = k + (k - 2) * h_dim1;
                            tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                    + (r__2 = h__[k + (k - 2) * h_dim1].imag, f2c_abs(r__2));
                        }
                        if(k >= *ktop + 3)
                        {
                            i__4 = k + (k - 3) * h_dim1;
                            tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                    + (r__2 = h__[k + (k - 3) * h_dim1].imag, f2c_abs(r__2));
                        }
                        if(k <= *kbot - 2)
                        {
                            i__4 = k + 2 + (k + 1) * h_dim1;
                            tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                    + (r__2 = h__[k + 2 + (k + 1) * h_dim1].imag, f2c_abs(r__2));
                        }
                        if(k <= *kbot - 3)
                        {
                            i__4 = k + 3 + (k + 1) * h_dim1;
                            tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                    + (r__2 = h__[k + 3 + (k + 1) * h_dim1].imag, f2c_abs(r__2));
                        }
                        if(k <= *kbot - 4)
                        {
                            i__4 = k + 4 + (k + 1) * h_dim1;
                            tst1 += (r__1 = h__[i__4].real, f2c_abs(r__1))
                                    + (r__2 = h__[k + 4 + (k + 1) * h_dim1].imag, f2c_abs(r__2));
                        }
                    }
                    i__4 = k + 1 + k * h_dim1;
                    i__5 = k + (k + 1) * h_dim1;
                    u1r = h__[i__4].real;
                    u1i = h__[i__4].imag;
                    u2r = h__[i__5].real;
                    u2i = h__[i__5].imag;

                    /* Computing MAX */
                    r__3 = smlnum;
                    r__4 = ulp * tst1; // , expr subst
                    if((r__1 = u1r, f2c_abs(r__1)) + (r__2 = u1i, f2c_abs(r__2))
                       <= fla_max(r__3, r__4))
                    {
                        /* Computing MAX */
                        r__5 = (r__1 = u1r, f2c_abs(r__1)) + (r__2 = u1i, f2c_abs(r__2));
                        r__6 = (r__3 = u2r, f2c_abs(r__3))
                               + (r__4 = u2i, f2c_abs(r__4)); // , expr subst
                        h12 = fla_max(r__5, r__6);
                        /* Computing MIN */
                        r__5 = (r__1 = u1r, f2c_abs(r__1)) + (r__2 = u1i, f2c_abs(r__2));
                        r__6 = (r__3 = u2r, f2c_abs(r__3))
                               + (r__4 = u2i, f2c_abs(r__4)); // , expr subst
                        h21 = fla_min(r__5, r__6);
                        i__4 = k + k * h_dim1;
                        i__5 = k + 1 + (k + 1) * h_dim1;
                        u1r = h__[i__4].real;
                        u1i = h__[i__4].imag;
                        u2r = h__[i__5].real;
                        u2i = h__[i__5].imag;
                        q__2.real = u1r - u2r;
                        q__2.imag = u1i - u2i; // , expr subst
                        q__1.real = q__2.real;
                        q__1.imag = q__2.imag; // , expr subst
                        /* Computing MAX */
                        r__5 = (r__1 = u2r, f2c_abs(r__1)) + (r__2 = u2i, f2c_abs(r__2));
                        r__6 = (r__3 = q__1.real, f2c_abs(r__3))
                               + (r__4 = q__1.imag, f2c_abs(r__4)); // , expr subst
                        h11 = fla_max(r__5, r__6);
                        q__2.real = u1r - u2r;
                        q__2.imag = u1i - u2i; // , expr subst
                        q__1.real = q__2.real;
                        q__1.imag = q__2.imag; // , expr subst
                        /* Computing MIN */
                        r__5 = (r__1 = u2r, f2c_abs(r__1)) + (r__2 = u2i, f2c_abs(r__2));
                        r__6 = (r__3 = q__1.real, f2c_abs(r__3))
                               + (r__4 = q__1.imag, f2c_abs(r__4)); // , expr subst
                        h22 = fla_min(r__5, r__6);
                        scl = h11 + h12;
                        tst2 = h22 * (h11 / scl);
                        /* Computing MAX */
                        r__1 = smlnum;
                        r__2 = ulp * tst2; // , expr subst
                        if(tst2 == 0.f || h21 * (h12 / scl) <= fla_max(r__1, r__2))
                        {
                            i__4 = k + 1 + k * h_dim1;
                            h__[i__4].real = 0.f;
                            h__[i__4].imag = 0.f; // , expr subst
                        }
                    }
                }
                /* L80: */
            }
            /* ==== Multiply H by reflections from the left ==== */
            if(accum)
            {
                jbot = fla_min(ndcol, *kbot);
            }
            else if(*wantt)
            {
                jbot = *n;
            }
            else
            {
                jbot = *kbot;
            }
            i__6 = mtop;
            for(m = mbot; m >= i__6; --m)
            {
                k = krcol + (m - 1 << 1);
                /* Computing MAX */
                i__4 = *ktop;
                i__5 = krcol + (m << 1); // , expr subst
                i__7 = jbot;

                v1r = v[m * v_dim1 + 1].real;
                v1i = v[m * v_dim1 + 1].imag;
                v2r = v[m * v_dim1 + 2].real;
                v2i = v[m * v_dim1 + 2].imag;
                v3r = v[m * v_dim1 + 3].real;
                v3i = v[m * v_dim1 + 3].imag;

                for(j = fla_max(i__4, i__5); j <= i__7; ++j)
                {
                    i__4 = k + 1 + j * h_dim1;
                    i__5 = k + 2 + j * h_dim1;
                    i__8 = k + 3 + j * h_dim1;
                    u1r = h__[i__4].real;
                    u1i = h__[i__4].imag;
                    u2r = h__[i__5].real;
                    u2i = h__[i__5].imag;
                    u3r = h__[i__8].real;
                    u3i = h__[i__8].imag;
                    q__5.real = v2r * u2r + v2i * u2i;
                    q__5.imag = v2r * u2i - v2i * u2r; // , expr subst
                    q__4.real = u1r + q__5.real;
                    q__4.imag = u1i + q__5.imag; // , expr subst
                    q__7.real = v3r * u3r + v3i * u3i;
                    q__7.imag = v3r * u3i - v3i * u3r; // , expr subst
                    q__3.real = q__4.real + q__7.real;
                    q__3.imag = q__4.imag + q__7.imag; // , expr subst
                    q__1.real = v1r * q__3.real + v1i * q__3.imag;
                    q__1.imag = v1r * q__3.imag - v1i * q__3.real; // , expr subst
                    refsum.real = q__1.real;
                    refsum.imag = q__1.imag; // , expr subst
                    q__1.real = u1r - refsum.real;
                    q__1.imag = u1i - refsum.imag; // , expr subst
                    h__[i__4].real = q__1.real;
                    h__[i__4].imag = q__1.imag; // , expr subst
                    q__2.real = refsum.real * v2r - refsum.imag * v2i;
                    q__2.imag = refsum.real * v2i + refsum.imag * v2r; // , expr subst
                    q__1.real = u2r - q__2.real;
                    q__1.imag = u2i - q__2.imag; // , expr subst
                    h__[i__5].real = q__1.real;
                    h__[i__5].imag = q__1.imag; // , expr subst
                    q__2.real = refsum.real * v3r - refsum.imag * v3i;
                    q__2.imag = refsum.real * v3i + refsum.imag * v3r; // , expr subst
                    q__1.real = u3r - q__2.real;
                    q__1.imag = u3i - q__2.imag; // , expr subst
                    h__[i__8].real = q__1.real;
                    h__[i__8].imag = q__1.imag; // , expr subst
                    /* L90: */
                }
                /* L100: */
            }
            /* ==== Accumulate orthogonal transformations. ==== */
            if(accum)
            {
                /* ==== Accumulate U. (If needed, update Z later */
                /* . with an efficient matrix-matrix */
                /* . multiply.) ==== */
                i__6 = mtop;
                for(m = mbot; m >= i__6; --m)
                {
                    k = krcol + (m - 1 << 1);
                    kms = k - incol;
                    /* Computing MAX */
                    i__7 = 1;
                    i__4 = *ktop - incol; // , expr subst
                    i2 = fla_max(i__7, i__4);
                    /* Computing MAX */
                    i__7 = i2;
                    i__4 = kms - (krcol - incol) + 1; // , expr subst
                    i2 = fla_max(i__7, i__4);
                    /* Computing MIN */
                    i__7 = kdu;
                    i__4 = krcol + (mbot - 1 << 1) - incol + 5; // , expr subst
                    i4 = fla_min(i__7, i__4);
                    i__7 = i4;

                    v1r = v[m * v_dim1 + 1].real;
                    v1i = v[m * v_dim1 + 1].imag;
                    v2r = v[m * v_dim1 + 2].real;
                    v2i = v[m * v_dim1 + 2].imag;
                    v3r = v[m * v_dim1 + 3].real;
                    v3i = v[m * v_dim1 + 3].imag;

                    for(j = i2; j <= i__7; ++j)
                    {
                        i__5 = j + (kms + 1) * u_dim1;
                        i__9 = j + (kms + 2) * u_dim1;
                        i__11 = j + (kms + 3) * u_dim1;
                        u1r = u[i__5].real;
                        u1i = u[i__5].imag;
                        u2r = u[i__9].real;
                        u2i = u[i__9].imag;
                        u3r = u[i__11].real;
                        u3i = u[i__11].imag;
                        q__4.real = v2r * u2r - v2i * u2i;
                        q__4.imag = v2r * u2i + v2i * u2r; // , expr subst
                        q__3.real = u1r + q__4.real;
                        q__3.imag = u1i + q__4.imag; // , expr subst
                        q__5.real = v3r * u3r - v3i * u3i;
                        q__5.imag = v3r * u3i + v3i * u3r; // , expr subst
                        q__2.real = q__3.real + q__5.real;
                        q__2.imag = q__3.imag + q__5.imag; // , expr subst
                        q__1.real = v1r * q__2.real - v1i * q__2.imag;
                        q__1.imag = v1r * q__2.imag + v1i * q__2.real; // , expr subst
                        refsum.real = q__1.real;
                        refsum.imag = q__1.imag; // , expr subst
                        q__1.real = u1r - refsum.real;
                        q__1.imag = u1i - refsum.imag; // , expr subst
                        u[i__5].real = q__1.real;
                        u[i__5].imag = q__1.imag; // , expr subst
                        q__2.real = refsum.real * v2r + refsum.imag * v2i;
                        q__2.imag = -refsum.real * v2i + refsum.imag * v2r; // , expr subst
                        q__1.real = u2r - q__2.real;
                        q__1.imag = u2i - q__2.imag; // , expr subst
                        u[i__9].real = q__1.real;
                        u[i__9].imag = q__1.imag; // , expr subst
                        q__2.real = refsum.real * v3r + refsum.imag * v3i;
                        q__2.imag = -refsum.real * v3i + refsum.imag * v3r; // , expr subst
                        q__1.real = u3r - q__2.real;
                        q__1.imag = u3i - q__2.imag; // , expr subst
                        u[i__11].real = q__1.real;
                        u[i__11].imag = q__1.imag; // , expr subst
                        /* L110: */
                    }
                    /* L120: */
                }
            }
            else if(*wantz)
            {
                /* ==== U is not accumulated, so update Z */
                /* . now by multiplying by reflections */
                /* . from the right. ==== */
                i__6 = mtop;
                for(m = mbot; m >= i__6; --m)
                {
                    k = krcol + (m - 1 << 1);
                    i__7 = *ihiz;
                    for(j = *iloz; j <= i__7; ++j)
                    {
                        i__4 = m * v_dim1 + 1;
                        i__5 = j + (k + 1) * z_dim1;
                        i__8 = m * v_dim1 + 2;
                        i__9 = j + (k + 2) * z_dim1;
                        q__4.real = v[i__8].real * z__[i__9].real - v[i__8].imag * z__[i__9].imag;
                        q__4.imag = v[i__8].real * z__[i__9].imag + v[i__8].imag * z__[i__9].real; // , expr subst
                        q__3.real = z__[i__5].real + q__4.real;
                        q__3.imag = z__[i__5].imag + q__4.imag; // , expr subst
                        i__10 = m * v_dim1 + 3;
                        i__11 = j + (k + 3) * z_dim1;
                        q__5.real = v[i__10].real * z__[i__11].real - v[i__10].imag * z__[i__11].imag;
                        q__5.imag
                            = v[i__10].real * z__[i__11].imag + v[i__10].imag * z__[i__11].real; // , expr subst
                        q__2.real = q__3.real + q__5.real;
                        q__2.imag = q__3.imag + q__5.imag; // , expr subst
                        q__1.real = v[i__4].real * q__2.real - v[i__4].imag * q__2.imag;
                        q__1.imag = v[i__4].real * q__2.imag + v[i__4].imag * q__2.real; // , expr subst
                        refsum.real = q__1.real;
                        refsum.imag = q__1.imag; // , expr subst
                        i__4 = j + (k + 1) * z_dim1;
                        i__5 = j + (k + 1) * z_dim1;
                        q__1.real = z__[i__5].real - refsum.real;
                        q__1.imag = z__[i__5].imag - refsum.imag; // , expr subst
                        z__[i__4].real = q__1.real;
                        z__[i__4].imag = q__1.imag; // , expr subst
                        i__4 = j + (k + 2) * z_dim1;
                        i__5 = j + (k + 2) * z_dim1;
                        q__3.real = v[m * v_dim1 + 2].real;
                        q__3.imag = -v[m * v_dim1 + 2].imag;
                        q__2.real = refsum.real * q__3.real - refsum.imag * q__3.imag;
                        q__2.imag = refsum.real * q__3.imag + refsum.imag * q__3.real; // , expr subst
                        q__1.real = z__[i__5].real - q__2.real;
                        q__1.imag = z__[i__5].imag - q__2.imag; // , expr subst
                        z__[i__4].real = q__1.real;
                        z__[i__4].imag = q__1.imag; // , expr subst
                        i__4 = j + (k + 3) * z_dim1;
                        i__5 = j + (k + 3) * z_dim1;
                        q__3.real = v[m * v_dim1 + 3].real;
                        q__3.imag = -v[m * v_dim1 + 3].imag;
                        q__2.real = refsum.real * q__3.real - refsum.imag * q__3.imag;
                        q__2.imag = refsum.real * q__3.imag + refsum.imag * q__3.real; // , expr subst
                        q__1.real = z__[i__5].real - q__2.real;
                        q__1.imag = z__[i__5].imag - q__2.imag; // , expr subst
                        z__[i__4].real = q__1.real;
                        z__[i__4].imag = q__1.imag; // , expr subst
                        /* L130: */
                    }
                    /* L140: */
                }
            }
            /* ==== End of near-the-diagonal bulge chase. ==== */
            /* L145: */
        }
        /* ==== Use U (if accumulated) to update far-from-diagonal */
        /* . entries in H. If required, use U to update Z as */
        /* . well. ==== */
        if(accum)
        {
            if(*wantt)
            {
                jtop = 1;
                jbot = *n;
            }
            else
            {
                jtop = *ktop;
                jbot = *kbot;
            }
            /* Computing MAX */
            i__3 = 1;
            i__6 = *ktop - incol; // , expr subst
            k1 = fla_max(i__3, i__6);
            /* Computing MAX */
            i__3 = 0;
            i__6 = ndcol - *kbot; // , expr subst
            nu = kdu - fla_max(i__3, i__6) - k1 + 1;
            /* ==== Horizontal Multiply ==== */
            i__3 = jbot;
            i__6 = *nh;
            for(jcol = fla_min(ndcol, *kbot) + 1; i__6 < 0 ? jcol >= i__3 : jcol <= i__3;
                jcol += i__6)
            {
                /* Computing MIN */
                i__7 = *nh;
                i__4 = jbot - jcol + 1; // , expr subst
                jlen = fla_min(i__7, i__4);
                aocl_blas_cgemm("C", "N", &nu, &jlen, &nu, &c_b2, &u[k1 + k1 * u_dim1], ldu,
                                &h__[incol + k1 + jcol * h_dim1], ldh, &c_b1, &wh[wh_offset], ldwh);
                aocl_lapack_clacpy("ALL", &nu, &jlen, &wh[wh_offset], ldwh,
                                   &h__[incol + k1 + jcol * h_dim1], ldh);
                /* L150: */
            }
            /* ==== Vertical multiply ==== */
            i__6 = fla_max(*ktop, incol) - 1;
            i__3 = *nv;
            for(jrow = jtop; i__3 < 0 ? jrow >= i__6 : jrow <= i__6; jrow += i__3)
            {
                /* Computing MIN */
                i__7 = *nv;
                i__4 = fla_max(*ktop, incol) - jrow; // , expr subst
                jlen = fla_min(i__7, i__4);
                aocl_blas_cgemm("N", "N", &jlen, &nu, &nu, &c_b2,
                                &h__[jrow + (incol + k1) * h_dim1], ldh, &u[k1 + k1 * u_dim1], ldu,
                                &c_b1, &wv[wv_offset], ldwv);
                aocl_lapack_clacpy("ALL", &jlen, &nu, &wv[wv_offset], ldwv,
                                   &h__[jrow + (incol + k1) * h_dim1], ldh);
                /* L160: */
            }
            /* ==== Z multiply (also vertical) ==== */
            if(*wantz)
            {
                i__3 = *ihiz;
                i__6 = *nv;
                for(jrow = *iloz; i__6 < 0 ? jrow >= i__3 : jrow <= i__3; jrow += i__6)
                {
                    /* Computing MIN */
                    i__7 = *nv;
                    i__4 = *ihiz - jrow + 1; // , expr subst
                    jlen = fla_min(i__7, i__4);
                    aocl_blas_cgemm("N", "N", &jlen, &nu, &nu, &c_b2,
                                    &z__[jrow + (incol + k1) * z_dim1], ldz, &u[k1 + k1 * u_dim1],
                                    ldu, &c_b1, &wv[wv_offset], ldwv);
                    aocl_lapack_clacpy("ALL", &jlen, &nu, &wv[wv_offset], ldwv,
                                       &z__[jrow + (incol + k1) * z_dim1], ldz);
                    /* L170: */
                }
            }
        }
        /* L180: */
    }
    /* ==== End of CLAQR5 ==== */
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* claqr5_ */

