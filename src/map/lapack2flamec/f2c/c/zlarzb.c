/* ../netlib/zlarzb.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
/* > \brief \b ZLARZB applies a block reflector or its conjugate-transpose to a general matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLARZB + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarzb.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarzb.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarzb.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, */
/* LDV, T, LDT, C, LDC, WORK, LDWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIRECT, SIDE, STOREV, TRANS */
/* INTEGER K, L, LDC, LDT, LDV, LDWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/* $ WORK( LDWORK, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARZB applies a complex block reflector H or its transpose H**H */
/* > to a complex distributed M-by-N C from the left or the right. */
/* > */
/* > Currently, only STOREV = 'R' and DIRECT = 'B' are supported. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply H or H**H from the Left */
/* > = 'R': apply H or H**H from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': apply H (No transpose) */
/* > = 'C': apply H**H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIRECT */
/* > \verbatim */
/* > DIRECT is CHARACTER*1 */
/* > Indicates how H is formed from a product of elementary */
/* > reflectors */
/* > = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet) */
/* > = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* > STOREV is CHARACTER*1 */
/* > Indicates how the vectors which define the elementary */
/* > reflectors are stored: */
/* > = 'C': Columnwise (not supported yet) */
/* > = 'R': Rowwise */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The order of the matrix T (= the number of elementary */
/* > reflectors whose product defines the block reflector). */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The number of columns of the matrix V containing the */
/* > meaningful part of the Householder reflectors. */
/* > If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is COMPLEX*16 array, dimension (LDV,NV). */
/* > If STOREV = 'C', NV = K;
if STOREV = 'R', NV = L. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. */
/* > If STOREV = 'C', LDV >= L;
if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is COMPLEX*16 array, dimension (LDT,K) */
/* > The triangular K-by-K matrix T in the representation of the */
/* > block reflector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (LDWORK,K) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* > LDWORK is INTEGER */
/* > The leading dimension of the array WORK. */
/* > If SIDE = 'L', LDWORK >= fla_max(1,N);
 */
/* > if SIDE = 'R', LDWORK >= fla_max(1,M). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void zlarzb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
             integer *k, integer *l, doublecomplex *v, integer *ldv, doublecomplex *t, integer *ldt,
             doublecomplex *c__, integer *ldc, doublecomplex *work, integer *ldwork)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlarzb inputs: side %c, trans %c, direct %c, storev %c, m %" FLA_IS
                      ", n %" FLA_IS ", k %" FLA_IS ", l %" FLA_IS ", ldv %" FLA_IS ", ldt %" FLA_IS
                      ", ldc %" FLA_IS ", ldwork % " FLA_IS "",
                      *side, *trans, *direct, *storev, *m, *n, *k, *l, *ldv, *ldt, *ldc, *ldwork);

    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, work_offset, i__1,
        i__2, i__3, i__4, i__5;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j, info;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *,
               integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *),
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        zlacgv_(integer *, doublecomplex *, integer *);
    char transt[1];
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    /* Function Body */
    if(*m <= 0 || *n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Check for currently supported options */
    info = 0;
    if(!lsame_(direct, "B", 1, 1))
    {
        info = -3;
    }
    else if(!lsame_(storev, "R", 1, 1))
    {
        info = -4;
    }
    if(info != 0)
    {
        i__1 = -info;
        xerbla_("ZLARZB", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(lsame_(trans, "N", 1, 1))
    {
        *(unsigned char *)transt = 'C';
    }
    else
    {
        *(unsigned char *)transt = 'N';
    }
    if(lsame_(side, "L", 1, 1))
    {
        /* Form H * C or H**H * C */
        /* W( 1:n, 1:k ) = C( 1:k, 1:n )**H */
        i__1 = *k;
        for(j = 1; j <= i__1; ++j)
        {
            zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
            /* L10: */
        }
        /* W( 1:n, 1:k ) = W( 1:n, 1:k ) + ... */
        /* C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T */
        if(*l > 0)
        {
            zgemm_("Transpose", "Conjugate transpose", n, k, l, &c_b1, &c__[*m - *l + 1 + c_dim1],
                   ldc, &v[v_offset], ldv, &c_b1, &work[work_offset], ldwork);
        }
        /* W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T or W( 1:m, 1:k ) * T */
        ztrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b1, &t[t_offset], ldt,
               &work[work_offset], ldwork);
        /* C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H */
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *k;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * c_dim1;
                i__4 = i__ + j * c_dim1;
                i__5 = j + i__ * work_dim1;
                z__1.r = c__[i__4].r - work[i__5].r;
                z__1.i = c__[i__4].i - work[i__5].i; // , expr subst
                c__[i__3].r = z__1.r;
                c__[i__3].i = z__1.i; // , expr subst
                /* L20: */
            }
            /* L30: */
        }
        /* C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ... */
        /* V( 1:k, 1:l )**H * W( 1:n, 1:k )**H */
        if(*l > 0)
        {
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemm_("Transpose", "Transpose", l, n, k, &z__1, &v[v_offset], ldv, &work[work_offset],
                   ldwork, &c_b1, &c__[*m - *l + 1 + c_dim1], ldc);
        }
    }
    else if(lsame_(side, "R", 1, 1))
    {
        /* Form C * H or C * H**H */
        /* W( 1:m, 1:k ) = C( 1:m, 1:k ) */
        i__1 = *k;
        for(j = 1; j <= i__1; ++j)
        {
            zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &c__1);
            /* L40: */
        }
        /* W( 1:m, 1:k ) = W( 1:m, 1:k ) + ... */
        /* C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H */
        if(*l > 0)
        {
            zgemm_("No transpose", "Transpose", m, k, l, &c_b1, &c__[(*n - *l + 1) * c_dim1 + 1],
                   ldc, &v[v_offset], ldv, &c_b1, &work[work_offset], ldwork);
        }
        /* W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T ) or */
        /* W( 1:m, 1:k ) * T**H */
        i__1 = *k;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *k - j + 1;
            zlacgv_(&i__2, &t[j + j * t_dim1], &c__1);
            /* L50: */
        }
        ztrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b1, &t[t_offset], ldt,
               &work[work_offset], ldwork);
        i__1 = *k;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *k - j + 1;
            zlacgv_(&i__2, &t[j + j * t_dim1], &c__1);
            /* L60: */
        }
        /* C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k ) */
        i__1 = *k;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * c_dim1;
                i__4 = i__ + j * c_dim1;
                i__5 = i__ + j * work_dim1;
                z__1.r = c__[i__4].r - work[i__5].r;
                z__1.i = c__[i__4].i - work[i__5].i; // , expr subst
                c__[i__3].r = z__1.r;
                c__[i__3].i = z__1.i; // , expr subst
                /* L70: */
            }
            /* L80: */
        }
        /* C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ... */
        /* W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) ) */
        i__1 = *l;
        for(j = 1; j <= i__1; ++j)
        {
            zlacgv_(k, &v[j * v_dim1 + 1], &c__1);
            /* L90: */
        }
        if(*l > 0)
        {
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemm_("No transpose", "No transpose", m, l, k, &z__1, &work[work_offset], ldwork,
                   &v[v_offset], ldv, &c_b1, &c__[(*n - *l + 1) * c_dim1 + 1], ldc);
        }
        i__1 = *l;
        for(j = 1; j <= i__1; ++j)
        {
            zlacgv_(k, &v[j * v_dim1 + 1], &c__1);
            /* L100: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLARZB */
}
/* zlarzb_ */
