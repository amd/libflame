/* ../netlib/clarft.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
static integer c__1 = 1;
/* > \brief \b CLARFT forms the triangular factor T of a block reflector H = I - vtvH */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARFT + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarft.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarft.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarft.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIRECT, STOREV */
/* INTEGER K, LDT, LDV, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX T( LDT, * ), TAU( * ), V( LDV, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFT forms the triangular factor T of a complex block reflector H */
/* > of order n, which is defined as a product of k elementary reflectors. */
/* > */
/* > If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
 */
/* > */
/* > If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. */
/* > */
/* > If STOREV = 'C', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th column of the array V, and */
/* > */
/* > H = I - V * T * V**H */
/* > */
/* > If STOREV = 'R', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th row of the array V, and */
/* > */
/* > H = I - V**H * T * V */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] DIRECT */
/* > \verbatim */
/* > DIRECT is CHARACTER*1 */
/* > Specifies the order in which the elementary reflectors are */
/* > multiplied to form the block reflector: */
/* > = 'F': H = H(1) H(2) . . . H(k) (Forward) */
/* > = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* > STOREV is CHARACTER*1 */
/* > Specifies how the vectors which define the elementary */
/* > reflectors are stored (see also Further Details): */
/* > = 'C': columnwise */
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
/* > \param[in] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension */
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
/* > TAU is COMPLEX array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX array, dimension (LDT,K) */
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
/* > \date December 2016 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The shape of the matrix V and the storage of the vectors which define */
/* > the H(i) is best illustrated by the following example with n = 5 and */
/* > k = 3. The elements equal to 1 are not stored. */
/* > */
/* > DIRECT = 'F' and STOREV = 'C': DIRECT = 'F' and STOREV = 'R': */
/* > */
/* > V = ( 1 ) V = ( 1 v1 v1 v1 v1 ) */
/* > ( v1 1 ) ( 1 v2 v2 v2 ) */
/* > ( v1 v2 1 ) ( 1 v3 v3 ) */
/* > ( v1 v2 v3 ) */
/* > ( v1 v2 v3 ) */
/* > */
/* > DIRECT = 'B' and STOREV = 'C': DIRECT = 'B' and STOREV = 'R': */
/* > */
/* > V = ( v1 v2 v3 ) V = ( v1 v1 1 ) */
/* > ( v1 v2 v3 ) ( v2 v2 v2 1 ) */
/* > ( 1 v2 v3 ) ( v3 v3 v3 v3 1 ) */
/* > ( 1 v3 ) */
/* > ( 1 ) */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void clarft_(char *direct, char *storev, integer *n, integer *k, complex *v, integer *ldv,
             complex *tau, complex *t, integer *ldt)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarft inputs: direct %c, storev %c, n %lld, k %lld, ldv %lld, ldt %lld",
             *direct, *storev, *n, *k, *ldv, *ldt);
#else
    snprintf(buffer, 256, "clarft inputs: direct %c, storev %c, n %d, k %d, ldv %d, ldt %d",
             *direct, *storev, *n, *k, *ldv, *ldt);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;
    /* Local variables */
    integer i__, j, prevlastv;
    extern /* Subroutine */
        void
        cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, complex *, integer *),
        cgemv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *,
               complex *, complex *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    integer lastv;
    extern /* Subroutine */
        void
        ctrmv_(char *, char *, char *, integer *, complex *, integer *, complex *, integer *),
        f90_exit_(void);
    /* -- LAPACK auxiliary routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
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
    /* Quick return if possible */
    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    /* Function Body */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(lsame_(direct, "F", 1, 1))
    {
        prevlastv = *n;
        i__1 = *k;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            prevlastv = fla_max(prevlastv, i__);
            i__2 = i__;
            if(tau[i__2].r == 0.f && tau[i__2].i == 0.f)
            {
                /* H(i) = I */
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = j + i__ * t_dim1;
                    t[i__3].r = 0.f;
                    t[i__3].i = 0.f; // , expr subst
                }
            }
            else
            {
                /* general case */
                if(lsame_(storev, "C", 1, 1))
                {
                    /* Skip any trailing zeros. */
                    i__2 = i__ + 1;
                    for(lastv = *n; lastv >= i__2; --lastv)
                    {
                        i__3 = lastv + i__ * v_dim1;
                        if(v[i__3].r != 0.f || v[i__3].i != 0.f)
                        {
                            break;
                        }
                    }
                    i__2 = i__ - 1;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = j + i__ * t_dim1;
                        i__4 = i__;
                        q__2.r = -tau[i__4].r;
                        q__2.i = -tau[i__4].i; // , expr subst
                        q__3.r = v[i__ + j * v_dim1].r;
                        q__3.i = -v[i__ + j * v_dim1].i;
                        q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                        q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                        t[i__3].r = q__1.r;
                        t[i__3].i = q__1.i; // , expr subst
                    }
                    j = fla_min(lastv, prevlastv);
                    /* T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i) */
                    i__2 = j - i__;
                    i__3 = i__ - 1;
                    i__4 = i__;
                    q__1.r = -tau[i__4].r;
                    q__1.i = -tau[i__4].i; // , expr subst
                    cgemv_("Conjugate transpose", &i__2, &i__3, &q__1, &v[i__ + 1 + v_dim1], ldv,
                           &v[i__ + 1 + i__ * v_dim1], &c__1, &c_b1, &t[i__ * t_dim1 + 1], &c__1);
                }
                else
                {
                    /* Skip any trailing zeros. */
                    i__2 = i__ + 1;
                    for(lastv = *n; lastv >= i__2; --lastv)
                    {
                        i__3 = i__ + lastv * v_dim1;
                        if(v[i__3].r != 0.f || v[i__3].i != 0.f)
                        {
                            break;
                        }
                    }
                    i__2 = i__ - 1;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = j + i__ * t_dim1;
                        i__4 = i__;
                        q__2.r = -tau[i__4].r;
                        q__2.i = -tau[i__4].i; // , expr subst
                        i__5 = j + i__ * v_dim1;
                        q__1.r = q__2.r * v[i__5].r - q__2.i * v[i__5].i;
                        q__1.i = q__2.r * v[i__5].i + q__2.i * v[i__5].r; // , expr subst
                        t[i__3].r = q__1.r;
                        t[i__3].i = q__1.i; // , expr subst
                    }
                    j = fla_min(lastv, prevlastv);
                    /* T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H */
                    i__2 = i__ - 1;
                    i__3 = j - i__;
                    i__4 = i__;
                    q__1.r = -tau[i__4].r;
                    q__1.i = -tau[i__4].i; // , expr subst
                    cgemm_("N", "C", &i__2, &c__1, &i__3, &q__1, &v[(i__ + 1) * v_dim1 + 1], ldv,
                           &v[i__ + (i__ + 1) * v_dim1], ldv, &c_b1, &t[i__ * t_dim1 + 1], ldt);
                }
                /* T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */
                i__2 = i__ - 1;
                ctrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt,
                       &t[i__ * t_dim1 + 1], &c__1);
                i__2 = i__ + i__ * t_dim1;
                i__3 = i__;
                t[i__2].r = tau[i__3].r;
                t[i__2].i = tau[i__3].i; // , expr subst
                if(i__ > 1)
                {
                    prevlastv = fla_max(prevlastv, lastv);
                }
                else
                {
                    prevlastv = lastv;
                }
            }
        }
    }
    else
    {
        prevlastv = 1;
        for(i__ = *k; i__ >= 1; --i__)
        {
            i__1 = i__;
            if(tau[i__1].r == 0.f && tau[i__1].i == 0.f)
            {
                /* H(i) = I */
                i__1 = *k;
                for(j = i__; j <= i__1; ++j)
                {
                    i__2 = j + i__ * t_dim1;
                    t[i__2].r = 0.f;
                    t[i__2].i = 0.f; // , expr subst
                }
            }
            else
            {
                /* general case */
                if(i__ < *k)
                {
                    if(lsame_(storev, "C", 1, 1))
                    {
                        /* Skip any leading zeros. */
                        i__1 = i__ - 1;
                        for(lastv = 1; lastv <= i__1; ++lastv)
                        {
                            i__2 = lastv + i__ * v_dim1;
                            if(v[i__2].r != 0.f || v[i__2].i != 0.f)
                            {
                                break;
                            }
                        }
                        i__1 = *k;
                        for(j = i__ + 1; j <= i__1; ++j)
                        {
                            i__2 = j + i__ * t_dim1;
                            i__3 = i__;
                            q__2.r = -tau[i__3].r;
                            q__2.i = -tau[i__3].i; // , expr subst
                            q__3.r = v[*n - *k + i__ + j * v_dim1].r;
                            q__3.i = -v[*n - *k + i__ + j * v_dim1].i;
                            q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                            q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                            t[i__2].r = q__1.r;
                            t[i__2].i = q__1.i; // , expr subst
                        }
                        j = fla_max(lastv, prevlastv);
                        /* T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i) */
                        i__1 = *n - *k + i__ - j;
                        i__2 = *k - i__;
                        i__3 = i__;
                        q__1.r = -tau[i__3].r;
                        q__1.i = -tau[i__3].i; // , expr subst
                        cgemv_("Conjugate transpose", &i__1, &i__2, &q__1,
                               &v[j + (i__ + 1) * v_dim1], ldv, &v[j + i__ * v_dim1], &c__1, &c_b1,
                               &t[i__ + 1 + i__ * t_dim1], &c__1);
                    }
                    else
                    {
                        /* Skip any leading zeros. */
                        i__1 = i__ - 1;
                        for(lastv = 1; lastv <= i__1; ++lastv)
                        {
                            i__2 = i__ + lastv * v_dim1;
                            if(v[i__2].r != 0.f || v[i__2].i != 0.f)
                            {
                                break;
                            }
                        }
                        i__1 = *k;
                        for(j = i__ + 1; j <= i__1; ++j)
                        {
                            i__2 = j + i__ * t_dim1;
                            i__3 = i__;
                            q__2.r = -tau[i__3].r;
                            q__2.i = -tau[i__3].i; // , expr subst
                            i__4 = j + (*n - *k + i__) * v_dim1;
                            q__1.r = q__2.r * v[i__4].r - q__2.i * v[i__4].i;
                            q__1.i = q__2.r * v[i__4].i + q__2.i * v[i__4].r; // , expr subst
                            t[i__2].r = q__1.r;
                            t[i__2].i = q__1.i; // , expr subst
                        }
                        j = fla_max(lastv, prevlastv);
                        /* T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H */
                        i__1 = *k - i__;
                        i__2 = *n - *k + i__ - j;
                        i__3 = i__;
                        q__1.r = -tau[i__3].r;
                        q__1.i = -tau[i__3].i; // , expr subst
                        cgemm_("N", "C", &i__1, &c__1, &i__2, &q__1, &v[i__ + 1 + j * v_dim1], ldv,
                               &v[i__ + j * v_dim1], ldv, &c_b1, &t[i__ + 1 + i__ * t_dim1], ldt);
                    }
                    /* T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i) */
                    i__1 = *k - i__;
                    ctrmv_("Lower", "No transpose", "Non-unit", &i__1,
                           &t[i__ + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ * t_dim1],
                           &c__1);
                    if(i__ > 1)
                    {
                        prevlastv = fla_min(prevlastv, lastv);
                    }
                    else
                    {
                        prevlastv = lastv;
                    }
                }
                i__1 = i__ + i__ * t_dim1;
                i__2 = i__;
                t[i__1].r = tau[i__2].r;
                t[i__1].i = tau[i__2].i; // , expr subst
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLARFT */
}
/* clarft_ */
