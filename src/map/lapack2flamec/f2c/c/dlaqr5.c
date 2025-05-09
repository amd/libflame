/* ./dlaqr5.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__3 = 3;
/* > \brief \b DLAQR5 performs a single small-bulge multi-shift QR sweep. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQR5 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr5.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr5.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr5.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, */
/* SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, */
/* LDU, NV, WV, LDWV, NH, WH, LDWH ) */
/* .. Scalar Arguments .. */
/* INTEGER IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, */
/* $ LDWH, LDWV, LDZ, N, NH, NSHFTS, NV */
/* LOGICAL WANTT, WANTZ */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION H( LDH, * ), SI( * ), SR( * ), U( LDU, * ), */
/* $ V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ), */
/* $ Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQR5, called by DLAQR0, performs a */
/* > single small-bulge multi-shift QR sweep. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTT */
/* > \verbatim */
/* > WANTT is LOGICAL */
/* > WANTT = .true. if the quasi-triangular Schur factor */
/* > is being computed. WANTT is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > WANTZ = .true. if the orthogonal Schur factor is being */
/* > computed. WANTZ is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] KACC22 */
/* > \verbatim */
/* > KACC22 is INTEGER with value 0, 1, or 2. */
/* > Specifies the computation mode of far-from-diagonal */
/* > orthogonal updates. */
/* > = 0: DLAQR5 does not accumulate reflections and does not */
/* > use matrix-matrix multiply to update far-from-diagonal */
/* > matrix entries. */
/* > = 1: DLAQR5 accumulates reflections and uses matrix-matrix */
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
/* > \param[in,out] SR */
/* > \verbatim */
/* > SR is DOUBLE PRECISION array, dimension (NSHFTS) */
/* > \endverbatim */
/* > */
/* > \param[in,out] SI */
/* > \verbatim */
/* > SI is DOUBLE PRECISION array, dimension (NSHFTS) */
/* > SR contains the real parts and SI contains the imaginary */
/* > parts of the NSHFTS shifts of origin that define the */
/* > multi-shift QR sweep. On output SR and SI may be */
/* > reordered. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* > H is DOUBLE PRECISION array, dimension (LDH,N) */
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
/* > Z is DOUBLE PRECISION array, dimension (LDZ,IHIZ) */
/* > If WANTZ = .TRUE., then the QR Sweep orthogonal */
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
/* > V is DOUBLE PRECISION array, dimension (LDV,NSHFTS/2) */
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
/* > U is DOUBLE PRECISION array, dimension (LDU,2*NSHFTS) */
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
/* > WV is DOUBLE PRECISION array, dimension (LDWV,2*NSHFTS) */
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
/* > WH is DOUBLE PRECISION array, dimension (LDWH,NH) */
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
/* > \ingroup laqr5 */
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
void dlaqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
             integer *kbot, integer *nshfts, doublereal *sr, doublereal *si, doublereal *h__,
             integer *ldh, integer *iloz, integer *ihiz, doublereal *z__, integer *ldz,
             doublereal *v, integer *ldv, doublereal *u, integer *ldu, integer *nv, doublereal *wv,
             integer *ldwv, integer *nh, doublereal *wh, integer *ldwh)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlaqr5 inputs: kacc22 %" FLA_IS ", n %" FLA_IS ", ktop %" FLA_IS
                      ", kbot %" FLA_IS ", nshfts %" FLA_IS ", ldh %" FLA_IS ", iloz %" FLA_IS
                      ", ihiz %" FLA_IS ", ldz %" FLA_IS ", ldv %" FLA_IS ", ldu %" FLA_IS
                      ", nv %" FLA_IS ", ldwv %" FLA_IS ", nh %" FLA_IS ", ldwh %" FLA_IS "",
                      *kacc22, *n, *ktop, *kbot, *nshfts, *ldh, *iloz, *ihiz, *ldz, *ldv, *ldu, *nv,
                      *ldwv, *nh, *ldwh);
    /* System generated locals */
    integer h_dim1, h_offset, u_dim1, u_offset, v_dim1, v_offset, wh_dim1, wh_offset, wv_dim1,
        wv_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3, d__4, d__5;
    /* Local variables */
    integer i__, j, k, m, i2, k1, i4;
    doublereal t1, t2, t3, h11, h12, h21, h22;
    integer m22, ns, nu;
    doublereal vt[3], scl;
    integer kdu, kms;
    doublereal ulp, tst1, tst2, beta;
    logical bmp22;
    integer jcol, jlen, jbot, mbot;
    doublereal swap;
    integer jtop, jrow, mtop;
    doublereal alpha;
    logical accum;
    extern /* Subroutine */
        void
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    integer ndcol, incol, krcol, nbmps;
    extern /* Subroutine */
        void
        dlaqr1_(integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
        void
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */
        void
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *);
    doublereal refsum, smlnum;
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
    /* .. Executable Statements .. */
    /* ==== If there are no shifts, then there is nothing to do. ==== */
    /* Parameter adjustments */
    --sr;
    --si;
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
    if(*nshfts < 2)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ==== If the active block is empty or 1-by-1, then there */
    /* . is nothing to do. ==== */
    if(*ktop >= *kbot)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ==== Shuffle shifts into pairs of real shifts and pairs */
    /* . of complex conjugate shifts assuming complex */
    /* . conjugate shifts are already adjacent to one */
    /* . another. ==== */
    i__1 = *nshfts - 2;
    for(i__ = 1; i__ <= i__1; i__ += 2)
    {
        if(si[i__] != -si[i__ + 1])
        {
            swap = sr[i__];
            sr[i__] = sr[i__ + 1];
            sr[i__ + 1] = sr[i__ + 2];
            sr[i__ + 2] = swap;
            swap = si[i__];
            si[i__] = si[i__ + 1];
            si[i__ + 1] = si[i__ + 2];
            si[i__ + 2] = swap;
        }
        /* L10: */
    }
    /* ==== NSHFTS is supposed to be even, but if it is odd, */
    /* . then simply reduce it by one. The shuffle above */
    /* . ensures that the dropped shift is real and that */
    /* . the remaining shifts are paired. ==== */
    ns = *nshfts - *nshfts % 2;
    /* ==== Machine constants for deflation ==== */
    safmin = dlamch_("SAFE MINIMUM");
    ulp = dlamch_("PRECISION");
    smlnum = safmin * ((doublereal)(*n) / ulp);
    /* ==== Use accumulated reflections to update far-from-diagonal */
    /* . entries ? ==== */
    accum = *kacc22 == 1 || *kacc22 == 2;
    /* ==== clear trash ==== */
    if(*ktop + 2 <= *kbot)
    {
        h__[*ktop + 2 + *ktop * h_dim1] = 0.;
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
            dlaset_("ALL", &kdu, &kdu, &c_b7, &c_b8, &u[u_offset], ldu);
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
                    dlaqr1_(&c__2, &h__[k + 1 + (k + 1) * h_dim1], ldh, &sr[(m22 << 1) - 1],
                            &si[(m22 << 1) - 1], &sr[m22 * 2], &si[m22 * 2], &v[m22 * v_dim1 + 1]);
                    beta = v[m22 * v_dim1 + 1];
                    dlarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 * v_dim1 + 1]);
                }
                else
                {
                    beta = h__[k + 1 + k * h_dim1];
                    v[m22 * v_dim1 + 2] = h__[k + 2 + k * h_dim1];
                    dlarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 * v_dim1 + 1]);
                    h__[k + 1 + k * h_dim1] = beta;
                    h__[k + 2 + k * h_dim1] = 0.;
                }
                /* ==== Perform update from right within */
                /* . computational window. ==== */
                t1 = v[m22 * v_dim1 + 1];
                t2 = t1 * v[m22 * v_dim1 + 2];
                /* Computing MIN */
                i__5 = *kbot;
                i__6 = k + 3; // , expr subst
                i__4 = fla_min(i__5, i__6);
                for(j = jtop; j <= i__4; ++j)
                {
                    refsum = h__[j + (k + 1) * h_dim1]
                             + v[m22 * v_dim1 + 2] * h__[j + (k + 2) * h_dim1];
                    h__[j + (k + 1) * h_dim1] -= refsum * t1;
                    h__[j + (k + 2) * h_dim1] -= refsum * t2;
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
                t1 = v[m22 * v_dim1 + 1];
                t2 = t1 * v[m22 * v_dim1 + 2];
                i__4 = jbot;
                for(j = k + 1; j <= i__4; ++j)
                {
                    refsum
                        = h__[k + 1 + j * h_dim1] + v[m22 * v_dim1 + 2] * h__[k + 2 + j * h_dim1];
                    h__[k + 1 + j * h_dim1] -= refsum * t1;
                    h__[k + 2 + j * h_dim1] -= refsum * t2;
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
                    if(h__[k + 1 + k * h_dim1] != 0.)
                    {
                        tst1 = (d__1 = h__[k + k * h_dim1], f2c_dabs(d__1))
                               + (d__2 = h__[k + 1 + (k + 1) * h_dim1], f2c_dabs(d__2));
                        if(tst1 == 0.)
                        {
                            if(k >= *ktop + 1)
                            {
                                tst1 += (d__1 = h__[k + (k - 1) * h_dim1], f2c_dabs(d__1));
                            }
                            if(k >= *ktop + 2)
                            {
                                tst1 += (d__1 = h__[k + (k - 2) * h_dim1], f2c_dabs(d__1));
                            }
                            if(k >= *ktop + 3)
                            {
                                tst1 += (d__1 = h__[k + (k - 3) * h_dim1], f2c_dabs(d__1));
                            }
                            if(k <= *kbot - 2)
                            {
                                tst1 += (d__1 = h__[k + 2 + (k + 1) * h_dim1], f2c_dabs(d__1));
                            }
                            if(k <= *kbot - 3)
                            {
                                tst1 += (d__1 = h__[k + 3 + (k + 1) * h_dim1], f2c_dabs(d__1));
                            }
                            if(k <= *kbot - 4)
                            {
                                tst1 += (d__1 = h__[k + 4 + (k + 1) * h_dim1], f2c_dabs(d__1));
                            }
                        }
                        /* Computing MAX */
                        d__2 = smlnum;
                        d__3 = ulp * tst1; // , expr subst
                        if((d__1 = h__[k + 1 + k * h_dim1], f2c_dabs(d__1)) <= fla_max(d__2, d__3))
                        {
                            /* Computing MAX */
                            d__3 = (d__1 = h__[k + 1 + k * h_dim1], f2c_dabs(d__1));
                            d__4 = (d__2 = h__[k + (k + 1) * h_dim1],
                                    f2c_dabs(d__2)); // , expr subst
                            h12 = fla_max(d__3, d__4);
                            /* Computing MIN */
                            d__3 = (d__1 = h__[k + 1 + k * h_dim1], f2c_dabs(d__1));
                            d__4 = (d__2 = h__[k + (k + 1) * h_dim1],
                                    f2c_dabs(d__2)); // , expr subst
                            h21 = fla_min(d__3, d__4);
                            /* Computing MAX */
                            d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], f2c_dabs(d__1));
                            d__4 = (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1],
                                    f2c_dabs(d__2)); // , expr subst
                            h11 = fla_max(d__3, d__4);
                            /* Computing MIN */
                            d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], f2c_dabs(d__1));
                            d__4 = (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1],
                                    f2c_dabs(d__2)); // , expr subst
                            h22 = fla_min(d__3, d__4);
                            scl = h11 + h12;
                            tst2 = h22 * (h11 / scl);
                            /* Computing MAX */
                            d__1 = smlnum;
                            d__2 = ulp * tst2; // , expr subst
                            if(tst2 == 0. || h21 * (h12 / scl) <= fla_max(d__1, d__2))
                            {
                                h__[k + 1 + k * h_dim1] = 0.;
                            }
                        }
                    }
                }
                /* ==== Accumulate orthogonal transformations. ==== */
                if(accum)
                {
                    kms = k - incol;
                    t1 = v[m22 * v_dim1 + 1];
                    t2 = t1 * v[m22 * v_dim1 + 2];
                    /* Computing MAX */
                    i__4 = 1;
                    i__5 = *ktop - incol; // , expr subst
                    i__6 = kdu;
                    for(j = fla_max(i__4, i__5); j <= i__6; ++j)
                    {
                        refsum = u[j + (kms + 1) * u_dim1]
                                 + v[m22 * v_dim1 + 2] * u[j + (kms + 2) * u_dim1];
                        u[j + (kms + 1) * u_dim1] -= refsum * t1;
                        u[j + (kms + 2) * u_dim1] -= refsum * t2;
                        /* L50: */
                    }
                }
                else if(*wantz)
                {
                    t1 = v[m22 * v_dim1 + 1];
                    t2 = t1 * v[m22 * v_dim1 + 2];
                    i__6 = *ihiz;
                    for(j = *iloz; j <= i__6; ++j)
                    {
                        refsum = z__[j + (k + 1) * z_dim1]
                                 + v[m22 * v_dim1 + 2] * z__[j + (k + 2) * z_dim1];
                        z__[j + (k + 1) * z_dim1] -= refsum * t1;
                        z__[j + (k + 2) * z_dim1] -= refsum * t2;
                        /* L60: */
                    }
                }
            }
            /* ==== Normal case: Chain of 3-by-3 reflections ==== */
            i__6 = mtop;
            for(m = mbot; m >= i__6; --m)
            {
                k = krcol + (m - 1 << 1);
                if(k == *ktop - 1)
                {
                    dlaqr1_(&c__3, &h__[*ktop + *ktop * h_dim1], ldh, &sr[(m << 1) - 1],
                            &si[(m << 1) - 1], &sr[m * 2], &si[m * 2], &v[m * v_dim1 + 1]);
                    alpha = v[m * v_dim1 + 1];
                    dlarfg_(&c__3, &alpha, &v[m * v_dim1 + 2], &c__1, &v[m * v_dim1 + 1]);
                }
                else
                {
                    /* ==== Perform delayed transformation of row below */
                    /* . Mth bulge. Exploit fact that first two elements */
                    /* . of row are actually zero. ==== */
                    t1 = v[m * v_dim1 + 1];
                    t2 = t1 * v[m * v_dim1 + 2];
                    t3 = t1 * v[m * v_dim1 + 3];
                    refsum = v[m * v_dim1 + 3] * h__[k + 3 + (k + 2) * h_dim1];
                    h__[k + 3 + k * h_dim1] = -refsum * t1;
                    h__[k + 3 + (k + 1) * h_dim1] = -refsum * t2;
                    h__[k + 3 + (k + 2) * h_dim1] -= refsum * t3;
                    /* ==== Calculate reflection to move */
                    /* . Mth bulge one step. ==== */
                    beta = h__[k + 1 + k * h_dim1];
                    v[m * v_dim1 + 2] = h__[k + 2 + k * h_dim1];
                    v[m * v_dim1 + 3] = h__[k + 3 + k * h_dim1];
                    dlarfg_(&c__3, &beta, &v[m * v_dim1 + 2], &c__1, &v[m * v_dim1 + 1]);
                    /* ==== A Bulge may collapse because of vigilant */
                    /* . deflation or destructive underflow. In the */
                    /* . underflow case, try the two-small-subdiagonals */
                    /* . trick to try to reinflate the bulge. ==== */
                    if(h__[k + 3 + k * h_dim1] != 0. || h__[k + 3 + (k + 1) * h_dim1] != 0.
                       || h__[k + 3 + (k + 2) * h_dim1] == 0.)
                    {
                        /* ==== Typical case: not collapsed (yet). ==== */
                        h__[k + 1 + k * h_dim1] = beta;
                        h__[k + 2 + k * h_dim1] = 0.;
                        h__[k + 3 + k * h_dim1] = 0.;
                    }
                    else
                    {
                        /* ==== Atypical case: collapsed. Attempt to */
                        /* . reintroduce ignoring H(K+1,K) and H(K+2,K). */
                        /* . If the fill resulting from the new */
                        /* . reflector is too large, then abandon it. */
                        /* . Otherwise, use the new one. ==== */
                        dlaqr1_(&c__3, &h__[k + 1 + (k + 1) * h_dim1], ldh, &sr[(m << 1) - 1],
                                &si[(m << 1) - 1], &sr[m * 2], &si[m * 2], vt);
                        alpha = vt[0];
                        dlarfg_(&c__3, &alpha, &vt[1], &c__1, vt);
                        t1 = vt[0];
                        t2 = t1 * vt[1];
                        t3 = t1 * vt[2];
                        refsum = h__[k + 1 + k * h_dim1] + vt[1] * h__[k + 2 + k * h_dim1];
                        if((d__1 = h__[k + 2 + k * h_dim1] - refsum * t2, f2c_dabs(d__1))
                               + (d__2 = refsum * t3, f2c_dabs(d__2))
                           > ulp
                                 * ((d__3 = h__[k + k * h_dim1], f2c_dabs(d__3))
                                    + (d__4 = h__[k + 1 + (k + 1) * h_dim1], f2c_dabs(d__4))
                                    + (d__5 = h__[k + 2 + (k + 2) * h_dim1], f2c_dabs(d__5))))
                        {
                            /* ==== Starting a new bulge here would */
                            /* . create non-negligible fill. Use */
                            /* . the old one with trepidation. ==== */
                            h__[k + 1 + k * h_dim1] = beta;
                            h__[k + 2 + k * h_dim1] = 0.;
                            h__[k + 3 + k * h_dim1] = 0.;
                        }
                        else
                        {
                            /* ==== Starting a new bulge here would */
                            /* . create only negligible fill. */
                            /* . Replace the old reflector with */
                            /* . the new one. ==== */
                            h__[k + 1 + k * h_dim1] -= refsum * t1;
                            h__[k + 2 + k * h_dim1] = 0.;
                            h__[k + 3 + k * h_dim1] = 0.;
                            v[m * v_dim1 + 1] = vt[0];
                            v[m * v_dim1 + 2] = vt[1];
                            v[m * v_dim1 + 3] = vt[2];
                        }
                    }
                }
                /* ==== Apply reflection from the right and */
                /* . the first column of update from the left. */
                /* . These updates are required for the vigilant */
                /* . deflation check. We still delay most of the */
                /* . updates from the left for efficiency. ==== */
                t1 = v[m * v_dim1 + 1];
                t2 = t1 * v[m * v_dim1 + 2];
                t3 = t1 * v[m * v_dim1 + 3];
                /* Computing MIN */
                i__5 = *kbot;
                i__7 = k + 3; // , expr subst
                i__4 = fla_min(i__5, i__7);
                for(j = jtop; j <= i__4; ++j)
                {
                    refsum = h__[j + (k + 1) * h_dim1]
                             + v[m * v_dim1 + 2] * h__[j + (k + 2) * h_dim1]
                             + v[m * v_dim1 + 3] * h__[j + (k + 3) * h_dim1];
                    h__[j + (k + 1) * h_dim1] -= refsum * t1;
                    h__[j + (k + 2) * h_dim1] -= refsum * t2;
                    h__[j + (k + 3) * h_dim1] -= refsum * t3;
                    /* L70: */
                }
                /* ==== Perform update from left for subsequent */
                /* . column. ==== */
                refsum = h__[k + 1 + (k + 1) * h_dim1]
                         + v[m * v_dim1 + 2] * h__[k + 2 + (k + 1) * h_dim1]
                         + v[m * v_dim1 + 3] * h__[k + 3 + (k + 1) * h_dim1];
                h__[k + 1 + (k + 1) * h_dim1] -= refsum * t1;
                h__[k + 2 + (k + 1) * h_dim1] -= refsum * t2;
                h__[k + 3 + (k + 1) * h_dim1] -= refsum * t3;
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
                if(h__[k + 1 + k * h_dim1] != 0.)
                {
                    tst1 = (d__1 = h__[k + k * h_dim1], f2c_dabs(d__1))
                           + (d__2 = h__[k + 1 + (k + 1) * h_dim1], f2c_dabs(d__2));
                    if(tst1 == 0.)
                    {
                        if(k >= *ktop + 1)
                        {
                            tst1 += (d__1 = h__[k + (k - 1) * h_dim1], f2c_dabs(d__1));
                        }
                        if(k >= *ktop + 2)
                        {
                            tst1 += (d__1 = h__[k + (k - 2) * h_dim1], f2c_dabs(d__1));
                        }
                        if(k >= *ktop + 3)
                        {
                            tst1 += (d__1 = h__[k + (k - 3) * h_dim1], f2c_dabs(d__1));
                        }
                        if(k <= *kbot - 2)
                        {
                            tst1 += (d__1 = h__[k + 2 + (k + 1) * h_dim1], f2c_dabs(d__1));
                        }
                        if(k <= *kbot - 3)
                        {
                            tst1 += (d__1 = h__[k + 3 + (k + 1) * h_dim1], f2c_dabs(d__1));
                        }
                        if(k <= *kbot - 4)
                        {
                            tst1 += (d__1 = h__[k + 4 + (k + 1) * h_dim1], f2c_dabs(d__1));
                        }
                    }
                    /* Computing MAX */
                    d__2 = smlnum;
                    d__3 = ulp * tst1; // , expr subst
                    if((d__1 = h__[k + 1 + k * h_dim1], f2c_dabs(d__1)) <= fla_max(d__2, d__3))
                    {
                        /* Computing MAX */
                        d__3 = (d__1 = h__[k + 1 + k * h_dim1], f2c_dabs(d__1));
                        d__4 = (d__2 = h__[k + (k + 1) * h_dim1], f2c_dabs(d__2)); // , expr subst
                        h12 = fla_max(d__3, d__4);
                        /* Computing MIN */
                        d__3 = (d__1 = h__[k + 1 + k * h_dim1], f2c_dabs(d__1));
                        d__4 = (d__2 = h__[k + (k + 1) * h_dim1], f2c_dabs(d__2)); // , expr subst
                        h21 = fla_min(d__3, d__4);
                        /* Computing MAX */
                        d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], f2c_dabs(d__1));
                        d__4 = (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1],
                                f2c_dabs(d__2)); // , expr subst
                        h11 = fla_max(d__3, d__4);
                        /* Computing MIN */
                        d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], f2c_dabs(d__1));
                        d__4 = (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1],
                                f2c_dabs(d__2)); // , expr subst
                        h22 = fla_min(d__3, d__4);
                        scl = h11 + h12;
                        tst2 = h22 * (h11 / scl);
                        /* Computing MAX */
                        d__1 = smlnum;
                        d__2 = ulp * tst2; // , expr subst
                        if(tst2 == 0. || h21 * (h12 / scl) <= fla_max(d__1, d__2))
                        {
                            h__[k + 1 + k * h_dim1] = 0.;
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
                t1 = v[m * v_dim1 + 1];
                t2 = t1 * v[m * v_dim1 + 2];
                t3 = t1 * v[m * v_dim1 + 3];
                /* Computing MAX */
                i__4 = *ktop;
                i__5 = krcol + (m << 1); // , expr subst
                i__7 = jbot;
                for(j = fla_max(i__4, i__5); j <= i__7; ++j)
                {
                    refsum = h__[k + 1 + j * h_dim1] + v[m * v_dim1 + 2] * h__[k + 2 + j * h_dim1]
                             + v[m * v_dim1 + 3] * h__[k + 3 + j * h_dim1];
                    h__[k + 1 + j * h_dim1] -= refsum * t1;
                    h__[k + 2 + j * h_dim1] -= refsum * t2;
                    h__[k + 3 + j * h_dim1] -= refsum * t3;
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
                    t1 = v[m * v_dim1 + 1];
                    t2 = t1 * v[m * v_dim1 + 2];
                    t3 = t1 * v[m * v_dim1 + 3];
                    i__7 = i4;
                    for(j = i2; j <= i__7; ++j)
                    {
                        refsum = u[j + (kms + 1) * u_dim1]
                                 + v[m * v_dim1 + 2] * u[j + (kms + 2) * u_dim1]
                                 + v[m * v_dim1 + 3] * u[j + (kms + 3) * u_dim1];
                        u[j + (kms + 1) * u_dim1] -= refsum * t1;
                        u[j + (kms + 2) * u_dim1] -= refsum * t2;
                        u[j + (kms + 3) * u_dim1] -= refsum * t3;
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
                    t1 = v[m * v_dim1 + 1];
                    t2 = t1 * v[m * v_dim1 + 2];
                    t3 = t1 * v[m * v_dim1 + 3];
                    i__7 = *ihiz;
                    for(j = *iloz; j <= i__7; ++j)
                    {
                        refsum = z__[j + (k + 1) * z_dim1]
                                 + v[m * v_dim1 + 2] * z__[j + (k + 2) * z_dim1]
                                 + v[m * v_dim1 + 3] * z__[j + (k + 3) * z_dim1];
                        z__[j + (k + 1) * z_dim1] -= refsum * t1;
                        z__[j + (k + 2) * z_dim1] -= refsum * t2;
                        z__[j + (k + 3) * z_dim1] -= refsum * t3;
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
                dgemm_("C", "N", &nu, &jlen, &nu, &c_b8, &u[k1 + k1 * u_dim1], ldu,
                       &h__[incol + k1 + jcol * h_dim1], ldh, &c_b7, &wh[wh_offset], ldwh);
                dlacpy_("ALL", &nu, &jlen, &wh[wh_offset], ldwh, &h__[incol + k1 + jcol * h_dim1],
                        ldh);
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
                dgemm_("N", "N", &jlen, &nu, &nu, &c_b8, &h__[jrow + (incol + k1) * h_dim1], ldh,
                       &u[k1 + k1 * u_dim1], ldu, &c_b7, &wv[wv_offset], ldwv);
                dlacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv, &h__[jrow + (incol + k1) * h_dim1],
                        ldh);
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
                    dgemm_("N", "N", &jlen, &nu, &nu, &c_b8, &z__[jrow + (incol + k1) * z_dim1],
                           ldz, &u[k1 + k1 * u_dim1], ldu, &c_b7, &wv[wv_offset], ldwv);
                    dlacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv,
                            &z__[jrow + (incol + k1) * z_dim1], ldz);
                    /* L170: */
                }
            }
        }
        /* L180: */
    }
    /* ==== End of DLAQR5 ==== */
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* dlaqr5_ */
