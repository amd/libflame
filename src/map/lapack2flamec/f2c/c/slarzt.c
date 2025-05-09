/* ../netlib/slarzt.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b8 = 0.f;
static integer c__1 = 1;
/* > \brief \b SLARZT forms the triangular factor T of a block reflector H = I - vtvH. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARZT + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarzt.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarzt.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarzt.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIRECT, STOREV */
/* INTEGER K, LDT, LDV, N */
/* .. */
/* .. Array Arguments .. */
/* REAL T( LDT, * ), TAU( * ), V( LDV, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARZT forms the triangular factor T of a real block reflector */
/* > H of order > n, which is defined as a product of k elementary */
/* > reflectors. */
/* > */
/* > If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
 */
/* > */
/* > If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. */
/* > */
/* > If STOREV = 'C', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th column of the array V, and */
/* > */
/* > H = I - V * T * V**T */
/* > */
/* > If STOREV = 'R', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th row of the array V, and */
/* > */
/* > H = I - V**T * T * V */
/* > */
/* > Currently, only STOREV = 'R' and DIRECT = 'B' are supported. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] DIRECT */
/* > \verbatim */
/* > DIRECT is CHARACTER*1 */
/* > Specifies the order in which the elementary reflectors are */
/* > multiplied to form the block reflector: */
/* > = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet) */
/* > = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* > STOREV is CHARACTER*1 */
/* > Specifies how the vectors which define the elementary */
/* > reflectors are stored (see also Further Details): */
/* > = 'C': columnwise (not supported yet) */
/* > = 'R': rowwise */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the block reflector H. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The order of the triangular factor T (= the number of */
/* > elementary reflectors). K >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* > V is REAL array, dimension */
/* > (LDV,K) if STOREV = 'C' */
/* > (LDV,N) if STOREV = 'R' */
/* > The matrix V. See further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. */
/* > If STOREV = 'C', LDV >= fla_max(1,N);
if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is REAL array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is REAL array, dimension (LDT,K) */
/* > The k by k triangular factor T of the block reflector. */
/* > If DIRECT = 'F', T is upper triangular;
if DIRECT = 'B', T is */
/* > lower triangular. The rest of the array is not used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The shape of the matrix V and the storage of the vectors which define */
/* > the H(i) is best illustrated by the following example with n = 5 and */
/* > k = 3. The elements equal to 1 are not stored;
the corresponding */
/* > array elements are modified but restored on exit. The rest of the */
/* > array is not used. */
/* > */
/* > DIRECT = 'F' and STOREV = 'C': DIRECT = 'F' and STOREV = 'R': */
/* > */
/* > ______V_____ */
/* > ( v1 v2 v3 ) / \ */
/* > ( v1 v2 v3 ) ( v1 v1 v1 v1 v1 . . . . 1 ) */
/* > V = ( v1 v2 v3 ) ( v2 v2 v2 v2 v2 . . . 1 ) */
/* > ( v1 v2 v3 ) ( v3 v3 v3 v3 v3 . . 1 ) */
/* > ( v1 v2 v3 ) */
/* > . . . */
/* > . . . */
/* > 1 . . */
/* > 1 . */
/* > 1 */
/* > */
/* > DIRECT = 'B' and STOREV = 'C': DIRECT = 'B' and STOREV = 'R': */
/* > */
/* > ______V_____ */
/* > 1 / \ */
/* > . 1 ( 1 . . . . v1 v1 v1 v1 v1 ) */
/* > . . 1 ( . 1 . . . v2 v2 v2 v2 v2 ) */
/* > . . . ( . . 1 . . v3 v3 v3 v3 v3 ) */
/* > . . . */
/* > ( v1 v2 v3 ) */
/* > ( v1 v2 v3 ) */
/* > V = ( v1 v2 v3 ) */
/* > ( v1 v2 v3 ) */
/* > ( v1 v2 v3 ) */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void slarzt_(char *direct, char *storev, integer *n, integer *k, real *v, integer *ldv, real *tau,
             real *t, integer *ldt)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slarzt inputs: direct %c, storev %c, n %" FLA_IS ", k %" FLA_IS
                      ", ldv %" FLA_IS ", ldt %" FLA_IS "",
                      *direct, *storev, *n, *k, *ldv, *ldt);
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1;
    real r__1;
    /* Local variables */
    integer i__, j, info;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *,
               real *, integer *),
        strmv_(char *, char *, char *, integer *, real *, integer *, real *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Check for currently supported options */
    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    /* Function Body */
    info = 0;
    if(!lsame_(direct, "B", 1, 1))
    {
        info = -1;
    }
    else if(!lsame_(storev, "R", 1, 1))
    {
        info = -2;
    }
    if(info != 0)
    {
        i__1 = -info;
        xerbla_("SLARZT", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    for(i__ = *k; i__ >= 1; --i__)
    {
        if(tau[i__] == 0.f)
        {
            /* H(i) = I */
            i__1 = *k;
            for(j = i__; j <= i__1; ++j)
            {
                t[j + i__ * t_dim1] = 0.f;
                /* L10: */
            }
        }
        else
        {
            /* general case */
            if(i__ < *k)
            {
                /* T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**T */
                i__1 = *k - i__;
                r__1 = -tau[i__];
                sgemv_("No transpose", &i__1, n, &r__1, &v[i__ + 1 + v_dim1], ldv, &v[i__ + v_dim1],
                       ldv, &c_b8, &t[i__ + 1 + i__ * t_dim1], &c__1);
                /* T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i) */
                i__1 = *k - i__;
                strmv_("Lower", "No transpose", "Non-unit", &i__1, &t[i__ + 1 + (i__ + 1) * t_dim1],
                       ldt, &t[i__ + 1 + i__ * t_dim1], &c__1);
            }
            t[i__ + i__ * t_dim1] = tau[i__];
        }
        /* L20: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLARZT */
}
/* slarzt_ */
