/* ../netlib/sorbdb.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SORBDB */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORBDB + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, */
/* X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, */
/* TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIGNS, TRANS */
/* INTEGER INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, */
/* $ Q */
/* .. */
/* .. Array Arguments .. */
/* REAL PHI( * ), THETA( * ) */
/* REAL TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), */
/* $ WORK( * ), X11( LDX11, * ), X12( LDX12, * ), */
/* $ X21( LDX21, * ), X22( LDX22, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORBDB simultaneously bidiagonalizes the blocks of an M-by-M */
/* > partitioned orthogonal matrix X: */
/* > */
/* > [ B11 | B12 0 0 ] */
/* > [ X11 | X12 ] [ P1 | ] [ 0 | 0 -I 0 ] [ Q1 | ]**T */
/* > X = [-----------] = [---------] [----------------] [---------] . */
/* > [ X21 | X22 ] [ | P2 ] [ B21 | B22 0 0 ] [ | Q2 ] */
/* > [ 0 | 0 0 I ] */
/* > */
/* > X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is */
/* > not the case, then X must be transposed and/or permuted. This can be */
/* > done in constant time using the TRANS and SIGNS options. See SORCSD */
/* > for details.) */
/* > */
/* > The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by- */
/* > (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are */
/* > represented implicitly by Householder vectors. */
/* > */
/* > B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented */
/* > implicitly by angles THETA, PHI. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > The number of columns in X11 and X21. 0 <= Q <= */
/* > MIN(P,M-P,M-Q). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* > X11 is REAL array, dimension (LDX11,Q) */
/* > On entry, the top-left block of the orthogonal matrix to be */
/* > reduced. On exit, the form depends on TRANS: */
/* > If TRANS = 'N', then */
/* > the columns of tril(X11) specify reflectors for P1, */
/* > the rows of triu(X11,1) specify reflectors for Q1;
 */
/* > else TRANS = 'T', and */
/* > the rows of triu(X11) specify reflectors for P1, */
/* > the columns of tril(X11,-1) specify reflectors for Q1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* > LDX11 is INTEGER */
/* > The leading dimension of X11. If TRANS = 'N', then LDX11 >= */
/* > P;
else LDX11 >= Q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X12 */
/* > \verbatim */
/* > X12 is REAL array, dimension (LDX12,M-Q) */
/* > On entry, the top-right block of the orthogonal matrix to */
/* > be reduced. On exit, the form depends on TRANS: */
/* > If TRANS = 'N', then */
/* > the rows of triu(X12) specify the first P reflectors for */
/* > Q2;
 */
/* > else TRANS = 'T', and */
/* > the columns of tril(X12) specify the first P reflectors */
/* > for Q2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX12 */
/* > \verbatim */
/* > LDX12 is INTEGER */
/* > The leading dimension of X12. If TRANS = 'N', then LDX12 >= */
/* > P;
else LDX11 >= M-Q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* > X21 is REAL array, dimension (LDX21,Q) */
/* > On entry, the bottom-left block of the orthogonal matrix to */
/* > be reduced. On exit, the form depends on TRANS: */
/* > If TRANS = 'N', then */
/* > the columns of tril(X21) specify reflectors for P2;
 */
/* > else TRANS = 'T', and */
/* > the rows of triu(X21) specify reflectors for P2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* > LDX21 is INTEGER */
/* > The leading dimension of X21. If TRANS = 'N', then LDX21 >= */
/* > M-P;
else LDX21 >= Q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X22 */
/* > \verbatim */
/* > X22 is REAL array, dimension (LDX22,M-Q) */
/* > On entry, the bottom-right block of the orthogonal matrix to */
/* > be reduced. On exit, the form depends on TRANS: */
/* > If TRANS = 'N', then */
/* > the rows of triu(X22(Q+1:M-P,P+1:M-Q)) specify the last */
/* > M-P-Q reflectors for Q2, */
/* > else TRANS = 'T', and */
/* > the columns of tril(X22(P+1:M-Q,Q+1:M-P)) specify the last */
/* > M-P-Q reflectors for P2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX22 */
/* > \verbatim */
/* > LDX22 is INTEGER */
/* > The leading dimension of X22. If TRANS = 'N', then LDX22 >= */
/* > M-P;
else LDX22 >= M-Q. */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* > THETA is REAL array, dimension (Q) */
/* > The entries of the bidiagonal blocks B11, B12, B21, B22 can */
/* > be computed from the angles THETA and PHI. See Further */
/* > Details. */
/* > \endverbatim */
/* > */
/* > \param[out] PHI */
/* > \verbatim */
/* > PHI is REAL array, dimension (Q-1) */
/* > The entries of the bidiagonal blocks B11, B12, B21, B22 can */
/* > be computed from the angles THETA and PHI. See Further */
/* > Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP1 */
/* > \verbatim */
/* > TAUP1 is REAL array, dimension (P) */
/* > The scalar factors of the elementary reflectors that define */
/* > P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* > TAUP2 is REAL array, dimension (M-P) */
/* > The scalar factors of the elementary reflectors that define */
/* > P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* > TAUQ1 is REAL array, dimension (Q) */
/* > The scalar factors of the elementary reflectors that define */
/* > Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ2 */
/* > \verbatim */
/* > TAUQ2 is REAL array, dimension (M-Q) */
/* > The scalar factors of the elementary reflectors that define */
/* > Q2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= M-Q. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup realOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The bidiagonal blocks B11, B12, B21, and B22 are represented */
/* > implicitly by angles THETA(1), ..., THETA(Q) and PHI(1), ..., */
/* > PHI(Q-1). B11 and B21 are upper bidiagonal, while B21 and B22 are */
/* > lower bidiagonal. Every entry in each bidiagonal band is a product */
/* > of a sine or cosine of a THETA with a sine or cosine of a PHI. See */
/* > [1] or SORCSD for details. */
/* > */
/* > P1, P2, Q1, and Q2 are represented as products of elementary */
/* > reflectors. See SORCSD for details on generating P1, P2, Q1, and Q2 */
/* > using SORGQR and SORGLQ. */
/* > \endverbatim */
/* > \par References: */
/* ================ */
/* > */
/* > [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* > Algorithms, 50(1):33-65, 2009. */
/* > */
/* ===================================================================== */
/* Subroutine */
void sorbdb_(char *trans, char *signs, integer *m, integer *p, integer *q, real *x11,
             integer *ldx11, real *x12, integer *ldx12, real *x21, integer *ldx21, real *x22,
             integer *ldx22, real *theta, real *phi, real *taup1, real *taup2, real *tauq1,
             real *tauq2, real *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorbdb inputs: trans %c, signs %c, m %" FLA_IS ", p %" FLA_IS ", q %" FLA_IS
                      ", ldx11 %" FLA_IS ", ldx12 %" FLA_IS ", ldx21 %" FLA_IS ", "
                      "ldx22 %" FLA_IS "",
                      *trans, *signs, *m, *p, *q, *ldx11, *ldx12, *ldx21, *ldx22);
    /* System generated locals */
    integer x11_dim1, x11_offset, x12_dim1, x12_offset, x21_dim1, x21_offset, x22_dim1, x22_offset,
        i__1, i__2, i__3;
    real r__1;
    /* Builtin functions */
    double cos(doublereal), sin(doublereal), atan2(doublereal, doublereal);
    /* Local variables */
    logical colmajor;
    integer lworkmin, lworkopt, i__;
    real z1, z2, z3, z4;
    extern real snrm2_(integer *, real *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *),
        slarf_(char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *),
        saxpy_(integer *, real *, real *, integer *, real *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical lquery;
    extern /* Subroutine */
        void
        slarfgp_(integer *, real *, real *, integer *, real *);
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ==================================================================== */
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
    --phi;
    --taup1;
    --taup2;
    --tauq1;
    --tauq2;
    --work;
    /* Function Body */
    *info = 0;
    colmajor = !lsame_(trans, "T", 1, 1);
    if(!lsame_(signs, "O", 1, 1))
    {
        z1 = 1.f;
        z2 = 1.f;
        z3 = 1.f;
        z4 = 1.f;
    }
    else
    {
        z1 = 1.f;
        z2 = -1.f;
        z3 = 1.f;
        z4 = -1.f;
    }
    lquery = *lwork == -1;
    if(*m < 0)
    {
        *info = -3;
    }
    else if(*p < 0 || *p > *m)
    {
        *info = -4;
    }
    else if(*q < 0 || *q > *p || *q > *m - *p || *q > *m - *q)
    {
        *info = -5;
    }
    else if(colmajor && *ldx11 < fla_max(1, *p))
    {
        *info = -7;
    }
    else if(!colmajor && *ldx11 < fla_max(1, *q))
    {
        *info = -7;
    }
    else if(colmajor && *ldx12 < fla_max(1, *p))
    {
        *info = -9;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        if(!colmajor && *ldx12 < fla_max(i__1, i__2))
        {
            *info = -9;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = *m - *p; // , expr subst
            if(colmajor && *ldx21 < fla_max(i__1, i__2))
            {
                *info = -11;
            }
            else if(!colmajor && *ldx21 < fla_max(1, *q))
            {
                *info = -11;
            }
            else /* if(complicated condition) */
            {
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m - *p; // , expr subst
                if(colmajor && *ldx22 < fla_max(i__1, i__2))
                {
                    *info = -13;
                }
                else /* if(complicated condition) */
                {
                    /* Computing MAX */
                    i__1 = 1;
                    i__2 = *m - *q; // , expr subst
                    if(!colmajor && *ldx22 < fla_max(i__1, i__2))
                    {
                        *info = -13;
                    }
                }
            }
        }
    }
    /* Compute workspace */
    if(*info == 0)
    {
        lworkopt = *m - *q;
        lworkmin = *m - *q;
        work[1] = (real)lworkopt;
        if(*lwork < lworkmin && !lquery)
        {
            *info = -21;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("xORBDB", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Handle column-major and row-major separately */
    if(colmajor)
    {
        /* Reduce columns 1, ..., Q of X11, X12, X21, and X22 */
        i__1 = *q;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if(i__ == 1)
            {
                i__2 = *p - i__ + 1;
                sscal_(&i__2, &z1, &x11[i__ + i__ * x11_dim1], &c__1);
            }
            else
            {
                i__2 = *p - i__ + 1;
                r__1 = z1 * cos(phi[i__ - 1]);
                sscal_(&i__2, &r__1, &x11[i__ + i__ * x11_dim1], &c__1);
                i__2 = *p - i__ + 1;
                r__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
                saxpy_(&i__2, &r__1, &x12[i__ + (i__ - 1) * x12_dim1], &c__1,
                       &x11[i__ + i__ * x11_dim1], &c__1);
            }
            if(i__ == 1)
            {
                i__2 = *m - *p - i__ + 1;
                sscal_(&i__2, &z2, &x21[i__ + i__ * x21_dim1], &c__1);
            }
            else
            {
                i__2 = *m - *p - i__ + 1;
                r__1 = z2 * cos(phi[i__ - 1]);
                sscal_(&i__2, &r__1, &x21[i__ + i__ * x21_dim1], &c__1);
                i__2 = *m - *p - i__ + 1;
                r__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
                saxpy_(&i__2, &r__1, &x22[i__ + (i__ - 1) * x22_dim1], &c__1,
                       &x21[i__ + i__ * x21_dim1], &c__1);
            }
            i__2 = *m - *p - i__ + 1;
            i__3 = *p - i__ + 1;
            theta[i__] = atan2(snrm2_(&i__2, &x21[i__ + i__ * x21_dim1], &c__1),
                               snrm2_(&i__3, &x11[i__ + i__ * x11_dim1], &c__1));
            if(*p > i__)
            {
                i__2 = *p - i__ + 1;
                slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + 1 + i__ * x11_dim1], &c__1,
                         &taup1[i__]);
            }
            else if(*p == i__)
            {
                i__2 = *p - i__ + 1;
                slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + i__ * x11_dim1], &c__1,
                         &taup1[i__]);
            }
            x11[i__ + i__ * x11_dim1] = 1.f;
            if(*m - *p > i__)
            {
                i__2 = *m - *p - i__ + 1;
                slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + i__ * x21_dim1], &c__1,
                         &taup2[i__]);
            }
            else if(*m - *p == i__)
            {
                i__2 = *m - *p - i__ + 1;
                slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * x21_dim1], &c__1,
                         &taup2[i__]);
            }
            x21[i__ + i__ * x21_dim1] = 1.f;
            if(*q > i__)
            {
                i__2 = *p - i__ + 1;
                i__3 = *q - i__;
                slarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &taup1[i__],
                       &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &work[1]);
            }
            if(*m - *q + 1 > i__)
            {
                i__2 = *p - i__ + 1;
                i__3 = *m - *q - i__ + 1;
                slarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &taup1[i__],
                       &x12[i__ + i__ * x12_dim1], ldx12, &work[1]);
            }
            if(*q > i__)
            {
                i__2 = *m - *p - i__ + 1;
                i__3 = *q - i__;
                slarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &taup2[i__],
                       &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &work[1]);
            }
            if(*m - *q + 1 > i__)
            {
                i__2 = *m - *p - i__ + 1;
                i__3 = *m - *q - i__ + 1;
                slarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &taup2[i__],
                       &x22[i__ + i__ * x22_dim1], ldx22, &work[1]);
            }
            if(i__ < *q)
            {
                i__2 = *q - i__;
                r__1 = -z1 * z3 * sin(theta[i__]);
                sscal_(&i__2, &r__1, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
                i__2 = *q - i__;
                r__1 = z2 * z3 * cos(theta[i__]);
                saxpy_(&i__2, &r__1, &x21[i__ + (i__ + 1) * x21_dim1], ldx21,
                       &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
            }
            i__2 = *m - *q - i__ + 1;
            r__1 = -z1 * z4 * sin(theta[i__]);
            sscal_(&i__2, &r__1, &x12[i__ + i__ * x12_dim1], ldx12);
            i__2 = *m - *q - i__ + 1;
            r__1 = z2 * z4 * cos(theta[i__]);
            saxpy_(&i__2, &r__1, &x22[i__ + i__ * x22_dim1], ldx22, &x12[i__ + i__ * x12_dim1],
                   ldx12);
            if(i__ < *q)
            {
                i__2 = *q - i__;
                i__3 = *m - *q - i__ + 1;
                phi[i__] = atan2(snrm2_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], ldx11),
                                 snrm2_(&i__3, &x12[i__ + i__ * x12_dim1], ldx12));
            }
            if(i__ < *q)
            {
                if(*q - i__ == 1)
                {
                    i__2 = *q - i__;
                    slarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1],
                             &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &tauq1[i__]);
                }
                else
                {
                    i__2 = *q - i__;
                    slarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1],
                             &x11[i__ + (i__ + 2) * x11_dim1], ldx11, &tauq1[i__]);
                }
                x11[i__ + (i__ + 1) * x11_dim1] = 1.f;
            }
            if(*q + i__ - 1 < *m)
            {
                if(*m - *q == i__)
                {
                    i__2 = *m - *q - i__ + 1;
                    slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * x12_dim1], ldx12,
                             &tauq2[i__]);
                }
                else
                {
                    i__2 = *m - *q - i__ + 1;
                    slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (i__ + 1) * x12_dim1],
                             ldx12, &tauq2[i__]);
                }
            }
            x12[i__ + i__ * x12_dim1] = 1.f;
            if(i__ < *q)
            {
                i__2 = *p - i__;
                i__3 = *q - i__;
                slarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &tauq1[i__],
                       &x11[i__ + 1 + (i__ + 1) * x11_dim1], ldx11, &work[1]);
                i__2 = *m - *p - i__;
                i__3 = *q - i__;
                slarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &tauq1[i__],
                       &x21[i__ + 1 + (i__ + 1) * x21_dim1], ldx21, &work[1]);
            }
            if(*p > i__)
            {
                i__2 = *p - i__;
                i__3 = *m - *q - i__ + 1;
                slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &tauq2[i__],
                       &x12[i__ + 1 + i__ * x12_dim1], ldx12, &work[1]);
            }
            if(*m - *p > i__)
            {
                i__2 = *m - *p - i__;
                i__3 = *m - *q - i__ + 1;
                slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &tauq2[i__],
                       &x22[i__ + 1 + i__ * x22_dim1], ldx22, &work[1]);
            }
        }
        /* Reduce columns Q + 1, ..., P of X12, X22 */
        i__1 = *p;
        for(i__ = *q + 1; i__ <= i__1; ++i__)
        {
            i__2 = *m - *q - i__ + 1;
            r__1 = -z1 * z4;
            sscal_(&i__2, &r__1, &x12[i__ + i__ * x12_dim1], ldx12);
            if(i__ >= *m - *q)
            {
                i__2 = *m - *q - i__ + 1;
                slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * x12_dim1], ldx12,
                         &tauq2[i__]);
            }
            else
            {
                i__2 = *m - *q - i__ + 1;
                slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (i__ + 1) * x12_dim1], ldx12,
                         &tauq2[i__]);
            }
            x12[i__ + i__ * x12_dim1] = 1.f;
            if(*p > i__)
            {
                i__2 = *p - i__;
                i__3 = *m - *q - i__ + 1;
                slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &tauq2[i__],
                       &x12[i__ + 1 + i__ * x12_dim1], ldx12, &work[1]);
            }
            if(*m - *p - *q >= 1)
            {
                i__2 = *m - *p - *q;
                i__3 = *m - *q - i__ + 1;
                slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &tauq2[i__],
                       &x22[*q + 1 + i__ * x22_dim1], ldx22, &work[1]);
            }
        }
        /* Reduce columns P + 1, ..., M - Q of X12, X22 */
        i__1 = *m - *p - *q;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = *m - *p - *q - i__ + 1;
            r__1 = z2 * z4;
            sscal_(&i__2, &r__1, &x22[*q + i__ + (*p + i__) * x22_dim1], ldx22);
            if(i__ == *m - *p - *q)
            {
                i__2 = *m - *p - *q - i__ + 1;
                slarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1],
                         &x22[*q + i__ + (*p + i__) * x22_dim1], ldx22, &tauq2[*p + i__]);
            }
            else
            {
                i__2 = *m - *p - *q - i__ + 1;
                slarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1],
                         &x22[*q + i__ + (*p + i__ + 1) * x22_dim1], ldx22, &tauq2[*p + i__]);
            }
            x22[*q + i__ + (*p + i__) * x22_dim1] = 1.f;
            if(i__ < *m - *p - *q)
            {
                i__2 = *m - *p - *q - i__;
                i__3 = *m - *p - *q - i__ + 1;
                slarf_("R", &i__2, &i__3, &x22[*q + i__ + (*p + i__) * x22_dim1], ldx22,
                       &tauq2[*p + i__], &x22[*q + i__ + 1 + (*p + i__) * x22_dim1], ldx22,
                       &work[1]);
            }
        }
    }
    else
    {
        /* Reduce columns 1, ..., Q of X11, X12, X21, X22 */
        i__1 = *q;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if(i__ == 1)
            {
                i__2 = *p - i__ + 1;
                sscal_(&i__2, &z1, &x11[i__ + i__ * x11_dim1], ldx11);
            }
            else
            {
                i__2 = *p - i__ + 1;
                r__1 = z1 * cos(phi[i__ - 1]);
                sscal_(&i__2, &r__1, &x11[i__ + i__ * x11_dim1], ldx11);
                i__2 = *p - i__ + 1;
                r__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
                saxpy_(&i__2, &r__1, &x12[i__ - 1 + i__ * x12_dim1], ldx12,
                       &x11[i__ + i__ * x11_dim1], ldx11);
            }
            if(i__ == 1)
            {
                i__2 = *m - *p - i__ + 1;
                sscal_(&i__2, &z2, &x21[i__ + i__ * x21_dim1], ldx21);
            }
            else
            {
                i__2 = *m - *p - i__ + 1;
                r__1 = z2 * cos(phi[i__ - 1]);
                sscal_(&i__2, &r__1, &x21[i__ + i__ * x21_dim1], ldx21);
                i__2 = *m - *p - i__ + 1;
                r__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
                saxpy_(&i__2, &r__1, &x22[i__ - 1 + i__ * x22_dim1], ldx22,
                       &x21[i__ + i__ * x21_dim1], ldx21);
            }
            i__2 = *m - *p - i__ + 1;
            i__3 = *p - i__ + 1;
            theta[i__] = atan2(snrm2_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21),
                               snrm2_(&i__3, &x11[i__ + i__ * x11_dim1], ldx11));
            i__2 = *p - i__ + 1;
            slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) * x11_dim1], ldx11,
                     &taup1[i__]);
            x11[i__ + i__ * x11_dim1] = 1.f;
            if(i__ == *m - *p)
            {
                i__2 = *m - *p - i__ + 1;
                slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * x21_dim1], ldx21,
                         &taup2[i__]);
            }
            else
            {
                i__2 = *m - *p - i__ + 1;
                slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 1) * x21_dim1], ldx21,
                         &taup2[i__]);
            }
            x21[i__ + i__ * x21_dim1] = 1.f;
            if(*q > i__)
            {
                i__2 = *q - i__;
                i__3 = *p - i__ + 1;
                slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &taup1[i__],
                       &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[1]);
            }
            if(*m - *q + 1 > i__)
            {
                i__2 = *m - *q - i__ + 1;
                i__3 = *p - i__ + 1;
                slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &taup1[i__],
                       &x12[i__ + i__ * x12_dim1], ldx12, &work[1]);
            }
            if(*q > i__)
            {
                i__2 = *q - i__;
                i__3 = *m - *p - i__ + 1;
                slarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &taup2[i__],
                       &x21[i__ + 1 + i__ * x21_dim1], ldx21, &work[1]);
            }
            if(*m - *q + 1 > i__)
            {
                i__2 = *m - *q - i__ + 1;
                i__3 = *m - *p - i__ + 1;
                slarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &taup2[i__],
                       &x22[i__ + i__ * x22_dim1], ldx22, &work[1]);
            }
            if(i__ < *q)
            {
                i__2 = *q - i__;
                r__1 = -z1 * z3 * sin(theta[i__]);
                sscal_(&i__2, &r__1, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
                i__2 = *q - i__;
                r__1 = z2 * z3 * cos(theta[i__]);
                saxpy_(&i__2, &r__1, &x21[i__ + 1 + i__ * x21_dim1], &c__1,
                       &x11[i__ + 1 + i__ * x11_dim1], &c__1);
            }
            i__2 = *m - *q - i__ + 1;
            r__1 = -z1 * z4 * sin(theta[i__]);
            sscal_(&i__2, &r__1, &x12[i__ + i__ * x12_dim1], &c__1);
            i__2 = *m - *q - i__ + 1;
            r__1 = z2 * z4 * cos(theta[i__]);
            saxpy_(&i__2, &r__1, &x22[i__ + i__ * x22_dim1], &c__1, &x12[i__ + i__ * x12_dim1],
                   &c__1);
            if(i__ < *q)
            {
                i__2 = *q - i__;
                i__3 = *m - *q - i__ + 1;
                phi[i__] = atan2(snrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &c__1),
                                 snrm2_(&i__3, &x12[i__ + i__ * x12_dim1], &c__1));
            }
            if(i__ < *q)
            {
                if(*q - i__ == 1)
                {
                    i__2 = *q - i__;
                    slarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ + 1 + i__ * x11_dim1],
                             &c__1, &tauq1[i__]);
                }
                else
                {
                    i__2 = *q - i__;
                    slarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ + 2 + i__ * x11_dim1],
                             &c__1, &tauq1[i__]);
                }
                x11[i__ + 1 + i__ * x11_dim1] = 1.f;
            }
            if(*m - *q > i__)
            {
                i__2 = *m - *q - i__ + 1;
                slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * x12_dim1], &c__1,
                         &tauq2[i__]);
            }
            else
            {
                i__2 = *m - *q - i__ + 1;
                slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * x12_dim1], &c__1,
                         &tauq2[i__]);
            }
            x12[i__ + i__ * x12_dim1] = 1.f;
            if(i__ < *q)
            {
                i__2 = *q - i__;
                i__3 = *p - i__;
                slarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &c__1, &tauq1[i__],
                       &x11[i__ + 1 + (i__ + 1) * x11_dim1], ldx11, &work[1]);
                i__2 = *q - i__;
                i__3 = *m - *p - i__;
                slarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &c__1, &tauq1[i__],
                       &x21[i__ + 1 + (i__ + 1) * x21_dim1], ldx21, &work[1]);
            }
            i__2 = *m - *q - i__ + 1;
            i__3 = *p - i__;
            slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &tauq2[i__],
                   &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &work[1]);
            if(*m - *p - i__ > 0)
            {
                i__2 = *m - *q - i__ + 1;
                i__3 = *m - *p - i__;
                slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &tauq2[i__],
                       &x22[i__ + (i__ + 1) * x22_dim1], ldx22, &work[1]);
            }
        }
        /* Reduce columns Q + 1, ..., P of X12, X22 */
        i__1 = *p;
        for(i__ = *q + 1; i__ <= i__1; ++i__)
        {
            i__2 = *m - *q - i__ + 1;
            r__1 = -z1 * z4;
            sscal_(&i__2, &r__1, &x12[i__ + i__ * x12_dim1], &c__1);
            i__2 = *m - *q - i__ + 1;
            slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * x12_dim1], &c__1,
                     &tauq2[i__]);
            x12[i__ + i__ * x12_dim1] = 1.f;
            if(*p > i__)
            {
                i__2 = *m - *q - i__ + 1;
                i__3 = *p - i__;
                slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &tauq2[i__],
                       &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &work[1]);
            }
            if(*m - *p - *q >= 1)
            {
                i__2 = *m - *q - i__ + 1;
                i__3 = *m - *p - *q;
                slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &tauq2[i__],
                       &x22[i__ + (*q + 1) * x22_dim1], ldx22, &work[1]);
            }
        }
        /* Reduce columns P + 1, ..., M - Q of X12, X22 */
        i__1 = *m - *p - *q;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = *m - *p - *q - i__ + 1;
            r__1 = z2 * z4;
            sscal_(&i__2, &r__1, &x22[*p + i__ + (*q + i__) * x22_dim1], &c__1);
            if(*m - *p - *q == i__)
            {
                i__2 = *m - *p - *q - i__ + 1;
                slarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1],
                         &x22[*p + i__ + (*q + i__) * x22_dim1], &c__1, &tauq2[*p + i__]);
                x22[*p + i__ + (*q + i__) * x22_dim1] = 1.f;
            }
            else
            {
                i__2 = *m - *p - *q - i__ + 1;
                slarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1],
                         &x22[*p + i__ + 1 + (*q + i__) * x22_dim1], &c__1, &tauq2[*p + i__]);
                x22[*p + i__ + (*q + i__) * x22_dim1] = 1.f;
                i__2 = *m - *p - *q - i__ + 1;
                i__3 = *m - *p - *q - i__;
                slarf_("L", &i__2, &i__3, &x22[*p + i__ + (*q + i__) * x22_dim1], &c__1,
                       &tauq2[*p + i__], &x22[*p + i__ + (*q + i__ + 1) * x22_dim1], ldx22,
                       &work[1]);
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SORBDB */
}
/* sorbdb_ */
