/* ../netlib/dorcsd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
static logical c_false = FALSE_;
/* > \brief \b DORCSD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORCSD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorcsd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorcsd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorcsd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, */
/* SIGNS, M, P, Q, X11, LDX11, X12, */
/* LDX12, X21, LDX21, X22, LDX22, THETA, */
/* U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, */
/* LDV2T, WORK, LWORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS */
/* INTEGER INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12, */
/* $ LDX21, LDX22, LWORK, M, P, Q */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION THETA( * ) */
/* DOUBLE PRECISION U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), */
/* $ V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ), */
/* $ X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, */
/* $ * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORCSD computes the CS decomposition of an M-by-M partitioned */
/* > orthogonal matrix X: */
/* > */
/* > [ I 0 0 | 0 0 0 ] */
/* > [ 0 C 0 | 0 -S 0 ] */
/* > [ X11 | X12 ] [ U1 | ] [ 0 0 0 | 0 0 -I ] [ V1 | ]**T */
/* > X = [-----------] = [---------] [---------------------] [---------] . */
/* > [ X21 | X22 ] [ | U2 ] [ 0 0 0 | I 0 0 ] [ | V2 ] */
/* > [ 0 S 0 | 0 C 0 ] */
/* > [ 0 0 I | 0 0 0 ] */
/* > */
/* > X11 is P-by-Q. The orthogonal matrices U1, U2, V1, and V2 are P-by-P, */
/* > (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are */
/* > R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in */
/* > which R = MIN(P,M-P,Q,M-Q). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU1 */
/* > \verbatim */
/* > JOBU1 is CHARACTER */
/* > = 'Y': U1 is computed;
 */
/* > otherwise: U1 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU2 */
/* > \verbatim */
/* > JOBU2 is CHARACTER */
/* > = 'Y': U2 is computed;
 */
/* > otherwise: U2 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV1T */
/* > \verbatim */
/* > JOBV1T is CHARACTER */
/* > = 'Y': V1T is computed;
 */
/* > otherwise: V1T is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV2T */
/* > \verbatim */
/* > JOBV2T is CHARACTER */
/* > = 'Y': V2T is computed;
 */
/* > otherwise: V2T is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER */
/* > = 'T': X, U1, U2, V1T, and V2T are stored in row-major */
/* > order;
 */
/* > otherwise: X, U1, U2, V1T, and V2T are stored in column- */
/* > major order. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGNS */
/* > \verbatim */
/* > SIGNS is CHARACTER */
/* > = 'O': The lower-left block is made nonpositive (the */
/* > "other" convention);
 */
/* > otherwise: The upper-right block is made nonpositive (the */
/* > "default" convention). */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows and columns in X. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows in X11 and X12. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* > Q is INTEGER */
/* > The number of columns in X11 and X21. 0 <= Q <= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* > X11 is DOUBLE PRECISION array, dimension (LDX11,Q) */
/* > On entry, part of the orthogonal matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* > LDX11 is INTEGER */
/* > The leading dimension of X11. LDX11 >= MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X12 */
/* > \verbatim */
/* > X12 is DOUBLE PRECISION array, dimension (LDX12,M-Q) */
/* > On entry, part of the orthogonal matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX12 */
/* > \verbatim */
/* > LDX12 is INTEGER */
/* > The leading dimension of X12. LDX12 >= MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* > X21 is DOUBLE PRECISION array, dimension (LDX21,Q) */
/* > On entry, part of the orthogonal matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* > LDX21 is INTEGER */
/* > The leading dimension of X11. LDX21 >= MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X22 */
/* > \verbatim */
/* > X22 is DOUBLE PRECISION array, dimension (LDX22,M-Q) */
/* > On entry, part of the orthogonal matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX22 */
/* > \verbatim */
/* > LDX22 is INTEGER */
/* > The leading dimension of X11. LDX22 >= MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* > THETA is DOUBLE PRECISION array, dimension (R), in which R = */
/* > MIN(P,M-P,Q,M-Q). */
/* > C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and */
/* > S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ). */
/* > \endverbatim */
/* > */
/* > \param[out] U1 */
/* > \verbatim */
/* > U1 is DOUBLE PRECISION array, dimension (P) */
/* > If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU1 */
/* > \verbatim */
/* > LDU1 is INTEGER */
/* > The leading dimension of U1. If JOBU1 = 'Y', LDU1 >= */
/* > MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[out] U2 */
/* > \verbatim */
/* > U2 is DOUBLE PRECISION array, dimension (M-P) */
/* > If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal */
/* > matrix U2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* > LDU2 is INTEGER */
/* > The leading dimension of U2. If JOBU2 = 'Y', LDU2 >= */
/* > MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] V1T */
/* > \verbatim */
/* > V1T is DOUBLE PRECISION array, dimension (Q) */
/* > If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal */
/* > matrix V1**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV1T */
/* > \verbatim */
/* > LDV1T is INTEGER */
/* > The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >= */
/* > MAX(1,Q). */
/* > \endverbatim */
/* > */
/* > \param[out] V2T */
/* > \verbatim */
/* > V2T is DOUBLE PRECISION array, dimension (M-Q) */
/* > If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) orthogonal */
/* > matrix V2**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV2T */
/* > \verbatim */
/* > LDV2T is INTEGER */
/* > The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >= */
/* > MAX(1,M-Q). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > If INFO > 0 on exit, WORK(2:R) contains the values PHI(1), */
/* > ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* > define the matrix in intermediate bidiagonal-block form */
/* > remaining after nonconvergence. INFO specifies the number */
/* > of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the work array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (M-MIN(P, M-P, Q, M-Q)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: DBBCSD did not converge. See the description of WORK */
/* > above for details. */
/* > \endverbatim */
/* > \par References: */
/* ================ */
/* > */
/* > [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* > Algorithms, 50(1):33-65, 2009. */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void dorcsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
             integer *m, integer *p, integer *q, doublereal *x11, integer *ldx11, doublereal *x12,
             integer *ldx12, doublereal *x21, integer *ldx21, doublereal *x22, integer *ldx22,
             doublereal *theta, doublereal *u1, integer *ldu1, doublereal *u2, integer *ldu2,
             doublereal *v1t, integer *ldv1t, doublereal *v2t, integer *ldv2t, doublereal *work,
             integer *lwork, integer *iwork, integer *info)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, v2t_dim1, v2t_offset,
        x11_dim1, x11_offset, x12_dim1, x12_offset, x21_dim1, x21_offset, x22_dim1, x22_offset,
        i__1, i__2, i__3, i__4, i__5, i__6;
    /* Local variables */
    logical colmajor;
    integer lworkmin, lworkopt, i__, j, childinfo, lbbcsdwork, lorbdbwork, lorglqwork, lorgqrwork,
        ib11d, ib11e, ib12d, ib12e, ib21d, ib21e, ib22d, ib22e, iphi;
    logical defaultsigns;
    extern logical lsame_(char *, char *, integer, integer);
    integer lbbcsdworkmin, itaup1, itaup2, itauq1, itauq2, lbbcsdworkopt;
    logical wantu1, wantu2;
    extern /* Subroutine */
        void
        dbbcsd_(char *, char *, char *, char *, char *, integer *, integer *, integer *,
                doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *, integer *);
    integer ibbcsd, lorbdbworkopt;
    extern /* Subroutine */
        void
        dorbdb_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *, integer *);
    integer iorbdb, lorglqworkmin, lorgqrworkmin;
    extern /* Subroutine */
        void
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *),
        dlapmr_(logical *, integer *, integer *, doublereal *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        dlapmt_(logical *, integer *, integer *, doublereal *, integer *, integer *);
    integer lorglqworkopt;
    extern /* Subroutine */
        void
        dorglq_fla(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                   doublereal *, integer *, integer *);
    integer lorgqrworkopt, iorglq;
    extern /* Subroutine */
        int
        lapack_dorgqr(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                      doublereal *, integer *, integer *);
    integer iorgqr;
    char signst[1], transt[1];
    logical lquery, wantv1t, wantv2t;
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* =================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input arguments */
    /* Parameter adjustments */
    x11_dim1 = *ldx11;
    x11_offset = 1 + x11_dim1;
    x11 -= x11_offset;
    x12_dim1 = *ldx12;
    x12_offset = 1 + x12_dim1;
    x12 -= x12_offset;
    x21_dim1 = *ldx21;
    x21_offset = 1 + x21_dim1;
    x21 -= x21_offset;
    x22_dim1 = *ldx22;
    x22_offset = 1 + x22_dim1;
    x22 -= x22_offset;
    --theta;
    u1_dim1 = *ldu1;
    u1_offset = 1 + u1_dim1;
    u1 -= u1_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    v1t_dim1 = *ldv1t;
    v1t_offset = 1 + v1t_dim1;
    v1t -= v1t_offset;
    v2t_dim1 = *ldv2t;
    v2t_offset = 1 + v2t_dim1;
    v2t -= v2t_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    wantu1 = lsame_(jobu1, "Y", 1, 1);
    wantu2 = lsame_(jobu2, "Y", 1, 1);
    wantv1t = lsame_(jobv1t, "Y", 1, 1);
    wantv2t = lsame_(jobv2t, "Y", 1, 1);
    colmajor = !lsame_(trans, "T", 1, 1);
    defaultsigns = !lsame_(signs, "O", 1, 1);
    lquery = *lwork == -1;
    iorgqr = 0;
    iorglq = 0;
    iorbdb = 0;
    ibbcsd = 0;
    ibbcsd = 0;
    itauq2 = 0;
    itauq1 = 0;
    itaup2 = 0;
    itaup1 = 0;
    ib22e = 0;
    ib22d = 0;
    ib21e = 0;
    ib21d = 0;
    ib12e = 0;
    ib12d = 0;
    ib11e = 0;
    ib11d = 0;
    if(*m < 0)
    {
        *info = -7;
    }
    else if(*p < 0 || *p > *m)
    {
        *info = -8;
    }
    else if(*q < 0 || *q > *m)
    {
        *info = -9;
    }
    else if(colmajor && *ldx11 < fla_max(1, *p))
    {
        *info = -11;
    }
    else if(!colmajor && *ldx11 < fla_max(1, *q))
    {
        *info = -11;
    }
    else if(colmajor && *ldx12 < fla_max(1, *p))
    {
        *info = -13;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        if(!colmajor && *ldx12 < fla_max(i__1, i__2))
        {
            *info = -13;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = *m - *p; // , expr subst
            if(colmajor && *ldx21 < fla_max(i__1, i__2))
            {
                *info = -15;
            }
            else if(!colmajor && *ldx21 < fla_max(1, *q))
            {
                *info = -15;
            }
            else /* if(complicated condition) */
            {
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m - *p; // , expr subst
                if(colmajor && *ldx22 < fla_max(i__1, i__2))
                {
                    *info = -17;
                }
                else /* if(complicated condition) */
                {
                    /* Computing MAX */
                    i__1 = 1;
                    i__2 = *m - *q; // , expr subst
                    if(!colmajor && *ldx22 < fla_max(i__1, i__2))
                    {
                        *info = -17;
                    }
                    else if(wantu1 && *ldu1 < *p)
                    {
                        *info = -20;
                    }
                    else if(wantu2 && *ldu2 < *m - *p)
                    {
                        *info = -22;
                    }
                    else if(wantv1t && *ldv1t < *q)
                    {
                        *info = -24;
                    }
                    else if(wantv2t && *ldv2t < *m - *q)
                    {
                        *info = -26;
                    }
                }
            }
        }
    }
    /* Work with transpose if convenient */
    /* Computing MIN */
    i__1 = *p;
    i__2 = *m - *p; // , expr subst
    /* Computing MIN */
    i__3 = *q;
    i__4 = *m - *q; // , expr subst
    if(*info == 0 && fla_min(i__1, i__2) < fla_min(i__3, i__4))
    {
        if(colmajor)
        {
            *(unsigned char *)transt = 'T';
        }
        else
        {
            *(unsigned char *)transt = 'N';
        }
        if(defaultsigns)
        {
            *(unsigned char *)signst = 'O';
        }
        else
        {
            *(unsigned char *)signst = 'D';
        }
        dorcsd_(jobv1t, jobv2t, jobu1, jobu2, transt, signst, m, q, p, &x11[x11_offset], ldx11,
                &x21[x21_offset], ldx21, &x12[x12_offset], ldx12, &x22[x22_offset], ldx22,
                &theta[1], &v1t[v1t_offset], ldv1t, &v2t[v2t_offset], ldv2t, &u1[u1_offset], ldu1,
                &u2[u2_offset], ldu2, &work[1], lwork, &iwork[1], info);
        return;
    }
    /* Work with permutation [ 0 I;
    I 0 ] * X * [ 0 I;
    I 0 ] if */
    /* convenient */
    if(*info == 0 && *m - *q < *q)
    {
        if(defaultsigns)
        {
            *(unsigned char *)signst = 'O';
        }
        else
        {
            *(unsigned char *)signst = 'D';
        }
        i__1 = *m - *p;
        i__2 = *m - *q;
        dorcsd_(jobu2, jobu1, jobv2t, jobv1t, trans, signst, m, &i__1, &i__2, &x22[x22_offset],
                ldx22, &x21[x21_offset], ldx21, &x12[x12_offset], ldx12, &x11[x11_offset], ldx11,
                &theta[1], &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &v2t[v2t_offset], ldv2t,
                &v1t[v1t_offset], ldv1t, &work[1], lwork, &iwork[1], info);
        return;
    }
    /* Compute workspace */
    if(*info == 0)
    {
        iphi = 2;
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        itaup1 = iphi + fla_max(i__1, i__2);
        itaup2 = itaup1 + fla_max(1, *p);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *p; // , expr subst
        itauq1 = itaup2 + fla_max(i__1, i__2);
        itauq2 = itauq1 + fla_max(1, *q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        iorgqr = itauq2 + fla_max(i__1, i__2);
        i__1 = *m - *q;
        i__2 = *m - *q;
        i__3 = *m - *q;
        /* Computing MAX */
        i__5 = 1;
        i__6 = *m - *q; // , expr subst
        i__4 = fla_max(i__5, i__6);
        lapack_dorgqr(&i__1, &i__2, &i__3, &u1[u1_offset], &i__4, &u1[u1_offset], &work[1], &c_n1,
                      &childinfo);
        lorgqrworkopt = (integer)work[1];
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        lorgqrworkmin = fla_max(i__1, i__2);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        iorglq = itauq2 + fla_max(i__1, i__2);
        i__1 = *m - *q;
        i__2 = *m - *q;
        i__3 = *m - *q;
        /* Computing MAX */
        i__5 = 1;
        i__6 = *m - *q; // , expr subst
        i__4 = fla_max(i__5, i__6);
        dorglq_fla(&i__1, &i__2, &i__3, &u1[u1_offset], &i__4, &u1[u1_offset], &work[1], &c_n1,
                   &childinfo);
        lorglqworkopt = (integer)work[1];
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        lorglqworkmin = fla_max(i__1, i__2);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        iorbdb = itauq2 + fla_max(i__1, i__2);
        dorbdb_(trans, signs, m, p, q, &x11[x11_offset], ldx11, &x12[x12_offset], ldx12,
                &x21[x21_offset], ldx21, &x22[x22_offset], ldx22, &theta[1], &v1t[v1t_offset],
                &u1[u1_offset], &u2[u2_offset], &v1t[v1t_offset], &v2t[v2t_offset], &work[1], &c_n1,
                &childinfo);
        lorbdbworkopt = (integer)work[1];
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        ib11d = itauq2 + fla_max(i__1, i__2);
        ib11e = ib11d + fla_max(1, *q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ib12d = ib11e + fla_max(i__1, i__2);
        ib12e = ib12d + fla_max(1, *q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ib21d = ib12e + fla_max(i__1, i__2);
        ib21e = ib21d + fla_max(1, *q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ib22d = ib21e + fla_max(i__1, i__2);
        ib22e = ib22d + fla_max(1, *q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ibbcsd = ib22e + fla_max(i__1, i__2);
        dbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, &theta[1], &theta[1], &u1[u1_offset],
                ldu1, &u2[u2_offset], ldu2, &v1t[v1t_offset], ldv1t, &v2t[v2t_offset], ldv2t,
                &u1[u1_offset], &u1[u1_offset], &u1[u1_offset], &u1[u1_offset], &u1[u1_offset],
                &u1[u1_offset], &u1[u1_offset], &u1[u1_offset], &work[1], &c_n1, &childinfo);
        lbbcsdworkopt = (integer)work[1];
        lbbcsdworkmin = lbbcsdworkopt;
        /* Computing MAX */
        i__1 = iorgqr + lorgqrworkopt, i__2 = iorglq + lorglqworkopt, i__1 = fla_max(i__1, i__2),
        i__2 = iorbdb + lorbdbworkopt;
        i__1 = fla_max(i__1, i__2);
        i__2 = ibbcsd + lbbcsdworkopt; // ; expr subst
        lworkopt = fla_max(i__1, i__2) - 1;
        /* Computing MAX */
        i__1 = iorgqr + lorgqrworkmin, i__2 = iorglq + lorglqworkmin, i__1 = fla_max(i__1, i__2),
        i__2 = iorbdb + lorbdbworkopt;
        i__1 = fla_max(i__1, i__2);
        i__2 = ibbcsd + lbbcsdworkmin; // ; expr subst
        lworkmin = fla_max(i__1, i__2) - 1;
        work[1] = (doublereal)fla_max(lworkopt, lworkmin);
        if(*lwork < lworkmin && !lquery)
        {
            *info = -22;
        }
        else
        {
            lorgqrwork = *lwork - iorgqr + 1;
            lorglqwork = *lwork - iorglq + 1;
            lorbdbwork = *lwork - iorbdb + 1;
            lbbcsdwork = *lwork - ibbcsd + 1;
        }
    }
    /* Abort if any illegal arguments */
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORCSD", &i__1, (ftnlen)6);
        return;
    }
    else if(lquery)
    {
        return;
    }
    /* Transform to bidiagonal block form */
    dorbdb_(trans, signs, m, p, q, &x11[x11_offset], ldx11, &x12[x12_offset], ldx12,
            &x21[x21_offset], ldx21, &x22[x22_offset], ldx22, &theta[1], &work[iphi], &work[itaup1],
            &work[itaup2], &work[itauq1], &work[itauq2], &work[iorbdb], &lorbdbwork, &childinfo);
    /* Accumulate Householder reflectors */
    if(colmajor)
    {
        if(wantu1 && *p > 0)
        {
            dlacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1);
            lapack_dorgqr(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[iorgqr], &lorgqrwork,
                          info);
        }
        if(wantu2 && *m - *p > 0)
        {
            i__1 = *m - *p;
            dlacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], ldu2);
            i__1 = *m - *p;
            i__2 = *m - *p;
            lapack_dorgqr(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &work[iorgqr],
                          &lorgqrwork, info);
        }
        if(wantv1t && *q > 0)
        {
            i__1 = *q - 1;
            i__2 = *q - 1;
            dlacpy_("U", &i__1, &i__2, &x11[(x11_dim1 << 1) + 1], ldx11, &v1t[(v1t_dim1 << 1) + 2],
                    ldv1t);
            v1t[v1t_dim1 + 1] = 1.;
            i__1 = *q;
            for(j = 2; j <= i__1; ++j)
            {
                v1t[j * v1t_dim1 + 1] = 0.;
                v1t[j + v1t_dim1] = 0.;
            }
            i__1 = *q - 1;
            i__2 = *q - 1;
            i__3 = *q - 1;
            dorglq_fla(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, &work[itauq1],
                       &work[iorglq], &lorglqwork, info);
        }
        if(wantv2t && *m - *q > 0)
        {
            i__1 = *m - *q;
            dlacpy_("U", p, &i__1, &x12[x12_offset], ldx12, &v2t[v2t_offset], ldv2t);
            if(*m - *p > *q)
            {
                i__1 = *m - *p - *q;
                i__2 = *m - *p - *q;
                dlacpy_("U", &i__1, &i__2, &x22[*q + 1 + (*p + 1) * x22_dim1], ldx22,
                        &v2t[*p + 1 + (*p + 1) * v2t_dim1], ldv2t);
            }
            if(*m > *q)
            {
                i__1 = *m - *q;
                i__2 = *m - *q;
                i__3 = *m - *q;
                dorglq_fla(&i__1, &i__2, &i__3, &v2t[v2t_offset], ldv2t, &work[itauq2],
                           &work[iorglq], &lorglqwork, info);
            }
        }
    }
    else
    {
        if(wantu1 && *p > 0)
        {
            dlacpy_("U", q, p, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1);
            dorglq_fla(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[iorglq], &lorglqwork,
                       info);
        }
        if(wantu2 && *m - *p > 0)
        {
            i__1 = *m - *p;
            dlacpy_("U", q, &i__1, &x21[x21_offset], ldx21, &u2[u2_offset], ldu2);
            i__1 = *m - *p;
            i__2 = *m - *p;
            dorglq_fla(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &work[iorglq],
                       &lorglqwork, info);
        }
        if(wantv1t && *q > 0)
        {
            i__1 = *q - 1;
            i__2 = *q - 1;
            dlacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &v1t[(v1t_dim1 << 1) + 2], ldv1t);
            v1t[v1t_dim1 + 1] = 1.;
            i__1 = *q;
            for(j = 2; j <= i__1; ++j)
            {
                v1t[j * v1t_dim1 + 1] = 0.;
                v1t[j + v1t_dim1] = 0.;
            }
            i__1 = *q - 1;
            i__2 = *q - 1;
            i__3 = *q - 1;
            lapack_dorgqr(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, &work[itauq1],
                          &work[iorgqr], &lorgqrwork, info);
        }
        if(wantv2t && *m - *q > 0)
        {
            i__1 = *m - *q;
            dlacpy_("L", &i__1, p, &x12[x12_offset], ldx12, &v2t[v2t_offset], ldv2t);
            i__1 = *m - *p - *q;
            i__2 = *m - *p - *q;
            dlacpy_("L", &i__1, &i__2, &x22[*p + 1 + (*q + 1) * x22_dim1], ldx22,
                    &v2t[*p + 1 + (*p + 1) * v2t_dim1], ldv2t);
            i__1 = *m - *q;
            i__2 = *m - *q;
            i__3 = *m - *q;
            lapack_dorgqr(&i__1, &i__2, &i__3, &v2t[v2t_offset], ldv2t, &work[itauq2],
                          &work[iorgqr], &lorgqrwork, info);
        }
    }
    /* Compute the CSD of the matrix in bidiagonal-block form */
    dbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, &theta[1], &work[iphi], &u1[u1_offset],
            ldu1, &u2[u2_offset], ldu2, &v1t[v1t_offset], ldv1t, &v2t[v2t_offset], ldv2t,
            &work[ib11d], &work[ib11e], &work[ib12d], &work[ib12e], &work[ib21d], &work[ib21e],
            &work[ib22d], &work[ib22e], &work[ibbcsd], &lbbcsdwork, info);
    /* Permute rows and columns to place identity submatrices in top- */
    /* left corner of (1,1)-block and/or bottom-right corner of (1,2)- */
    /* block and/or bottom-right corner of (2,1)-block and/or top-left */
    /* corner of (2,2)-block */
    if(*q > 0 && wantu2)
    {
        i__1 = *q;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            iwork[i__] = *m - *p - *q + i__;
        }
        i__1 = *m - *p;
        for(i__ = *q + 1; i__ <= i__1; ++i__)
        {
            iwork[i__] = i__ - *q;
        }
        if(colmajor)
        {
            i__1 = *m - *p;
            i__2 = *m - *p;
            dlapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
        }
        else
        {
            i__1 = *m - *p;
            i__2 = *m - *p;
            dlapmr_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
        }
    }
    if(*m > 0 && wantv2t)
    {
        i__1 = *p;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            iwork[i__] = *m - *p - *q + i__;
        }
        i__1 = *m - *q;
        for(i__ = *p + 1; i__ <= i__1; ++i__)
        {
            iwork[i__] = i__ - *p;
        }
        if(!colmajor)
        {
            i__1 = *m - *q;
            i__2 = *m - *q;
            dlapmt_(&c_false, &i__1, &i__2, &v2t[v2t_offset], ldv2t, &iwork[1]);
        }
        else
        {
            i__1 = *m - *q;
            i__2 = *m - *q;
            dlapmr_(&c_false, &i__1, &i__2, &v2t[v2t_offset], ldv2t, &iwork[1]);
        }
    }
    return;
    /* End DORCSD */
}
/* dorcsd_ */
