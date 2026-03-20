/* ../netlib/claqr1.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H a nd specified shifts. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQR1 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr1.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr1.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr1.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQR1( N, H, LDH, S1, S2, V ) */
/* .. Scalar Arguments .. */
/* COMPLEX S1, S2 */
/* INTEGER LDH, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX H( LDH, * ), V( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given a 2-by-2 or 3-by-3 matrix H, CLAQR1 sets v to a */
/* > scalar multiple of the first column of the product */
/* > */
/* > (*) K = (H - s1*I)*(H - s2*I) */
/* > */
/* > scaling to avoid overflows and most underflows. */
/* > */
/* > This is useful for starting double implicit shift bulges */
/* > in the QR algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > Order of the matrix H. N must be either 2 or 3. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* > H is COMPLEX array, dimension (LDH,N) */
/* > The 2-by-2 or 3-by-3 matrix H in (*). */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > The leading dimension of H as declared in */
/* > the calling procedure. LDH >= N */
/* > \endverbatim */
/* > */
/* > \param[in] S1 */
/* > \verbatim */
/* > S1 is COMPLEX */
/* > \endverbatim */
/* > */
/* > \param[in] S2 */
/* > \verbatim */
/* > S2 is COMPLEX */
/* > */
/* > S1 and S2 are the shifts defining K in (*) above. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension (N) */
/* > A scalar multiple of the first column of the */
/* > matrix K in (*). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Karen Braman and Ralph Byers, Department of Mathematics, */
/* > University of Kansas, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void claqr1_(aocl_int_t *n, scomplex *h__, aocl_int_t *ldh, scomplex *s1, scomplex *s2, scomplex *v)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_claqr1(n, h__, ldh, s1, s2, v);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldh_64 = *ldh;

    aocl_lapack_claqr1(&n_64, h__, &ldh_64, s1, s2, v);
#endif
}

void aocl_lapack_claqr1(aocl_int64_t *n, scomplex *h__, aocl_int64_t *ldh, scomplex *s1, scomplex *s2,
                        scomplex *v)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "claqr1 inputs: n %lld, ldh %lld", *n, *ldh);
#else
    snprintf(buffer, 256, "claqr1 inputs: n %d, ldh %d", *n, *ldh);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t h_dim1, h_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    scomplex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8;
    /* Local variables */
    real s;
    scomplex h21s, h31s;
    /* -- LAPACK auxiliary routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ================================================================ */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --v;
    /* Function Body */
    if(*n != 2 && *n != 3)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*n == 2)
    {
        i__1 = h_dim1 + 1;
        q__2.real = h__[i__1].real - s2->real;
        q__2.imag = h__[i__1].imag - s2->imag; // , expr subst
        q__1.real = q__2.real;
        q__1.imag = q__2.imag; // , expr subst
        i__2 = h_dim1 + 2;
        s = (r__1 = q__1.real, f2c_abs(r__1)) + (r__2 = q__1.imag, f2c_abs(r__2))
            + ((r__3 = h__[i__2].real, f2c_abs(r__3)) + (r__4 = h__[i__2].imag, f2c_abs(r__4)));
        if(s == 0.f)
        {
            v[1].real = 0.f;
            v[1].imag = 0.f; // , expr subst
            v[2].real = 0.f;
            v[2].imag = 0.f; // , expr subst
        }
        else
        {
            q__1.real = h__[i__2].real / s;
            q__1.imag = h__[i__2].imag / s; // , expr subst
            h21s.real = q__1.real;
            h21s.imag = q__1.imag; // , expr subst
            i__2 = (h_dim1 << 1) + 1;
            q__2.real = h21s.real * h__[i__2].real - h21s.imag * h__[i__2].imag;
            q__2.imag = h21s.real * h__[i__2].imag + h21s.imag * h__[i__2].real; // , expr subst
            q__4.real = h__[i__1].real - s1->real;
            q__4.imag = h__[i__1].imag - s1->imag; // , expr subst
            q__6.real = h__[i__1].real - s2->real;
            q__6.imag = h__[i__1].imag - s2->imag; // , expr subst
            q__5.real = q__6.real / s;
            q__5.imag = q__6.imag / s; // , expr subst
            q__3.real = q__4.real * q__5.real - q__4.imag * q__5.imag;
            q__3.imag = q__4.real * q__5.imag + q__4.imag * q__5.real; // , expr subst
            q__1.real = q__2.real + q__3.real;
            q__1.imag = q__2.imag + q__3.imag; // , expr subst
            v[1].real = q__1.real;
            v[1].imag = q__1.imag; // , expr subst
            i__2 = (h_dim1 << 1) + 2;
            q__4.real = h__[i__1].real + h__[i__2].real;
            q__4.imag = h__[i__1].imag + h__[i__2].imag; // , expr subst
            q__3.real = q__4.real - s1->real;
            q__3.imag = q__4.imag - s1->imag; // , expr subst
            q__2.real = q__3.real - s2->real;
            q__2.imag = q__3.imag - s2->imag; // , expr subst
            q__1.real = h21s.real * q__2.real - h21s.imag * q__2.imag;
            q__1.imag = h21s.real * q__2.imag + h21s.imag * q__2.real; // , expr subst
            v[2].real = q__1.real;
            v[2].imag = q__1.imag; // , expr subst
        }
    }
    else
    {
        i__1 = h_dim1 + 1;
        q__2.real = h__[i__1].real - s2->real;
        q__2.imag = h__[i__1].imag - s2->imag; // , expr subst
        q__1.real = q__2.real;
        q__1.imag = q__2.imag; // , expr subst
        i__2 = h_dim1 + 2;
        i__3 = h_dim1 + 3;
        s = (r__1 = q__1.real, f2c_abs(r__1)) + (r__2 = q__1.imag, f2c_abs(r__2))
            + ((r__3 = h__[i__2].real, f2c_abs(r__3)) + (r__4 = h__[i__2].imag, f2c_abs(r__4)))
            + ((r__5 = h__[i__3].real, f2c_abs(r__5)) + (r__6 = h__[i__3].imag, f2c_abs(r__6)));
        if(s == 0.f)
        {
            v[1].real = 0.f;
            v[1].imag = 0.f; // , expr subst
            v[2].real = 0.f;
            v[2].imag = 0.f; // , expr subst
            v[3].real = 0.f;
            v[3].imag = 0.f; // , expr subst
        }
        else
        {
            q__1.real = h__[i__2].real / s;
            q__1.imag = h__[i__2].imag / s; // , expr subst
            h21s.real = q__1.real;
            h21s.imag = q__1.imag; // , expr subst
            q__1.real = h__[i__3].real / s;
            q__1.imag = h__[i__3].imag / s; // , expr subst
            h31s.real = q__1.real;
            h31s.imag = q__1.imag; // , expr subst
            q__4.real = h__[i__1].real - s1->real;
            q__4.imag = h__[i__1].imag - s1->imag; // , expr subst
            q__6.real = h__[i__1].real - s2->real;
            q__6.imag = h__[i__1].imag - s2->imag; // , expr subst
            q__5.real = q__6.real / s;
            q__5.imag = q__6.imag / s; // , expr subst
            q__3.real = q__4.real * q__5.real - q__4.imag * q__5.imag;
            q__3.imag = q__4.real * q__5.imag + q__4.imag * q__5.real; // , expr subst
            i__3 = (h_dim1 << 1) + 1;
            q__7.real = h__[i__3].real * h21s.real - h__[i__3].imag * h21s.imag;
            q__7.imag = h__[i__3].real * h21s.imag + h__[i__3].imag * h21s.real; // , expr subst
            q__2.real = q__3.real + q__7.real;
            q__2.imag = q__3.imag + q__7.imag; // , expr subst
            i__4 = h_dim1 * 3 + 1;
            q__8.real = h__[i__4].real * h31s.real - h__[i__4].imag * h31s.imag;
            q__8.imag = h__[i__4].real * h31s.imag + h__[i__4].imag * h31s.real; // , expr subst
            q__1.real = q__2.real + q__8.real;
            q__1.imag = q__2.imag + q__8.imag; // , expr subst
            v[1].real = q__1.real;
            v[1].imag = q__1.imag; // , expr subst
            i__2 = (h_dim1 << 1) + 2;
            q__5.real = h__[i__1].real + h__[i__2].real;
            q__5.imag = h__[i__1].imag + h__[i__2].imag; // , expr subst
            q__4.real = q__5.real - s1->real;
            q__4.imag = q__5.imag - s1->imag; // , expr subst
            q__3.real = q__4.real - s2->real;
            q__3.imag = q__4.imag - s2->imag; // , expr subst
            q__2.real = h21s.real * q__3.real - h21s.imag * q__3.imag;
            q__2.imag = h21s.real * q__3.imag + h21s.imag * q__3.real; // , expr subst
            i__3 = h_dim1 * 3 + 2;
            q__6.real = h__[i__3].real * h31s.real - h__[i__3].imag * h31s.imag;
            q__6.imag = h__[i__3].real * h31s.imag + h__[i__3].imag * h31s.real; // , expr subst
            q__1.real = q__2.real + q__6.real;
            q__1.imag = q__2.imag + q__6.imag; // , expr subst
            v[2].real = q__1.real;
            v[2].imag = q__1.imag; // , expr subst
            i__1 = h_dim1 + 1;
            i__2 = h_dim1 * 3 + 3;
            q__5.real = h__[i__1].real + h__[i__2].real;
            q__5.imag = h__[i__1].imag + h__[i__2].imag; // , expr subst
            q__4.real = q__5.real - s1->real;
            q__4.imag = q__5.imag - s1->imag; // , expr subst
            q__3.real = q__4.real - s2->real;
            q__3.imag = q__4.imag - s2->imag; // , expr subst
            q__2.real = h31s.real * q__3.real - h31s.imag * q__3.imag;
            q__2.imag = h31s.real * q__3.imag + h31s.imag * q__3.real; // , expr subst
            i__3 = (h_dim1 << 1) + 3;
            q__6.real = h21s.real * h__[i__3].real - h21s.imag * h__[i__3].imag;
            q__6.imag = h21s.real * h__[i__3].imag + h21s.imag * h__[i__3].real; // , expr subst
            q__1.real = q__2.real + q__6.real;
            q__1.imag = q__2.imag + q__6.imag; // , expr subst
            v[3].real = q__1.real;
            v[3].imag = q__1.imag; // , expr subst
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* claqr1_ */
