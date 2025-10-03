/* ../netlib/zlaqr1.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H a nd specified shifts. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAQR1 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr1.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr1.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr1.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V ) */
/* .. Scalar Arguments .. */
/* COMPLEX*16 S1, S2 */
/* INTEGER LDH, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 H( LDH, * ), V( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given a 2-by-2 or 3-by-3 matrix H, ZLAQR1 sets v to a */
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
/* > H is COMPLEX*16 array, dimension (LDH,N) */
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
/* > S1 is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[in] S2 */
/* > \verbatim */
/* > S2 is COMPLEX*16 */
/* > */
/* > S1 and S2 are the shifts defining K in (*) above. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX*16 array, dimension (N) */
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
/* > \ingroup complex16OTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Karen Braman and Ralph Byers, Department of Mathematics, */
/* > University of Kansas, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlaqr1_(aocl_int_t *n, dcomplex *h__, aocl_int_t *ldh, dcomplex *s1,
             dcomplex *s2, dcomplex *v)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlaqr1(n, h__, ldh, s1, s2, v);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldh_64 = *ldh;

    aocl_lapack_zlaqr1(&n_64, h__, &ldh_64, s1, s2, v);
#endif
}

void aocl_lapack_zlaqr1(aocl_int64_t *n, dcomplex *h__, aocl_int64_t *ldh, dcomplex *s1,
                        dcomplex *s2, dcomplex *v)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlaqr1 inputs: n %" FLA_IS ", ldh %" FLA_IS "", *n, *ldh);
    /* System generated locals */
    aocl_int64_t h_dim1, h_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    dcomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;
    /* Builtin functions */
    double d_imag(dcomplex *);
    /* Local variables */
    doublereal s;
    dcomplex h21s, h31s;
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
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 2)
    {
        i__1 = h_dim1 + 1;
        z__2.real = h__[i__1].real - s2->real;
        z__2.imag = h__[i__1].imag - s2->imag; // , expr subst
        z__1.real = z__2.real;
        z__1.imag = z__2.imag; // , expr subst
        i__2 = h_dim1 + 2;
        s = (d__1 = z__1.real, f2c_abs(d__1)) + (d__2 = d_imag(&z__1), f2c_abs(d__2))
            + ((d__3 = h__[i__2].real, f2c_abs(d__3))
               + (d__4 = d_imag(&h__[h_dim1 + 2]), f2c_abs(d__4)));
        if(s == 0.)
        {
            v[1].real = 0.;
            v[1].imag = 0.; // , expr subst
            v[2].real = 0.;
            v[2].imag = 0.; // , expr subst
        }
        else
        {
            i__1 = h_dim1 + 2;
            z__1.real = h__[i__1].real / s;
            z__1.imag = h__[i__1].imag / s; // , expr subst
            h21s.real = z__1.real;
            h21s.imag = z__1.imag; // , expr subst
            i__1 = (h_dim1 << 1) + 1;
            z__2.real = h21s.real * h__[i__1].real - h21s.imag * h__[i__1].imag;
            z__2.imag = h21s.real * h__[i__1].imag + h21s.imag * h__[i__1].real; // , expr subst
            i__2 = h_dim1 + 1;
            z__4.real = h__[i__2].real - s1->real;
            z__4.imag = h__[i__2].imag - s1->imag; // , expr subst
            i__3 = h_dim1 + 1;
            z__6.real = h__[i__3].real - s2->real;
            z__6.imag = h__[i__3].imag - s2->imag; // , expr subst
            z__5.real = z__6.real / s;
            z__5.imag = z__6.imag / s; // , expr subst
            z__3.real = z__4.real * z__5.real - z__4.imag * z__5.imag;
            z__3.imag = z__4.real * z__5.imag + z__4.imag * z__5.real; // , expr subst
            z__1.real = z__2.real + z__3.real;
            z__1.imag = z__2.imag + z__3.imag; // , expr subst
            v[1].real = z__1.real;
            v[1].imag = z__1.imag; // , expr subst
            i__1 = h_dim1 + 1;
            i__2 = (h_dim1 << 1) + 2;
            z__4.real = h__[i__1].real + h__[i__2].real;
            z__4.imag = h__[i__1].imag + h__[i__2].imag; // , expr subst
            z__3.real = z__4.real - s1->real;
            z__3.imag = z__4.imag - s1->imag; // , expr subst
            z__2.real = z__3.real - s2->real;
            z__2.imag = z__3.imag - s2->imag; // , expr subst
            z__1.real = h21s.real * z__2.real - h21s.imag * z__2.imag;
            z__1.imag = h21s.real * z__2.imag + h21s.imag * z__2.real; // , expr subst
            v[2].real = z__1.real;
            v[2].imag = z__1.imag; // , expr subst
        }
    }
    else
    {
        i__1 = h_dim1 + 1;
        z__2.real = h__[i__1].real - s2->real;
        z__2.imag = h__[i__1].imag - s2->imag; // , expr subst
        z__1.real = z__2.real;
        z__1.imag = z__2.imag; // , expr subst
        i__2 = h_dim1 + 2;
        i__3 = h_dim1 + 3;
        s = (d__1 = z__1.real, f2c_abs(d__1)) + (d__2 = d_imag(&z__1), f2c_abs(d__2))
            + ((d__3 = h__[i__2].real, f2c_abs(d__3))
               + (d__4 = d_imag(&h__[h_dim1 + 2]), f2c_abs(d__4)))
            + ((d__5 = h__[i__3].real, f2c_abs(d__5))
               + (d__6 = d_imag(&h__[h_dim1 + 3]), f2c_abs(d__6)));
        if(s == 0.)
        {
            v[1].real = 0.;
            v[1].imag = 0.; // , expr subst
            v[2].real = 0.;
            v[2].imag = 0.; // , expr subst
            v[3].real = 0.;
            v[3].imag = 0.; // , expr subst
        }
        else
        {
            i__1 = h_dim1 + 2;
            z__1.real = h__[i__1].real / s;
            z__1.imag = h__[i__1].imag / s; // , expr subst
            h21s.real = z__1.real;
            h21s.imag = z__1.imag; // , expr subst
            i__1 = h_dim1 + 3;
            z__1.real = h__[i__1].real / s;
            z__1.imag = h__[i__1].imag / s; // , expr subst
            h31s.real = z__1.real;
            h31s.imag = z__1.imag; // , expr subst
            i__1 = h_dim1 + 1;
            z__4.real = h__[i__1].real - s1->real;
            z__4.imag = h__[i__1].imag - s1->imag; // , expr subst
            i__2 = h_dim1 + 1;
            z__6.real = h__[i__2].real - s2->real;
            z__6.imag = h__[i__2].imag - s2->imag; // , expr subst
            z__5.real = z__6.real / s;
            z__5.imag = z__6.imag / s; // , expr subst
            z__3.real = z__4.real * z__5.real - z__4.imag * z__5.imag;
            z__3.imag = z__4.real * z__5.imag + z__4.imag * z__5.real; // , expr subst
            i__3 = (h_dim1 << 1) + 1;
            z__7.real = h__[i__3].real * h21s.real - h__[i__3].imag * h21s.imag;
            z__7.imag = h__[i__3].real * h21s.imag + h__[i__3].imag * h21s.real; // , expr subst
            z__2.real = z__3.real + z__7.real;
            z__2.imag = z__3.imag + z__7.imag; // , expr subst
            i__4 = h_dim1 * 3 + 1;
            z__8.real = h__[i__4].real * h31s.real - h__[i__4].imag * h31s.imag;
            z__8.imag = h__[i__4].real * h31s.imag + h__[i__4].imag * h31s.real; // , expr subst
            z__1.real = z__2.real + z__8.real;
            z__1.imag = z__2.imag + z__8.imag; // , expr subst
            v[1].real = z__1.real;
            v[1].imag = z__1.imag; // , expr subst
            i__1 = h_dim1 + 1;
            i__2 = (h_dim1 << 1) + 2;
            z__5.real = h__[i__1].real + h__[i__2].real;
            z__5.imag = h__[i__1].imag + h__[i__2].imag; // , expr subst
            z__4.real = z__5.real - s1->real;
            z__4.imag = z__5.imag - s1->imag; // , expr subst
            z__3.real = z__4.real - s2->real;
            z__3.imag = z__4.imag - s2->imag; // , expr subst
            z__2.real = h21s.real * z__3.real - h21s.imag * z__3.imag;
            z__2.imag = h21s.real * z__3.imag + h21s.imag * z__3.real; // , expr subst
            i__3 = h_dim1 * 3 + 2;
            z__6.real = h__[i__3].real * h31s.real - h__[i__3].imag * h31s.imag;
            z__6.imag = h__[i__3].real * h31s.imag + h__[i__3].imag * h31s.real; // , expr subst
            z__1.real = z__2.real + z__6.real;
            z__1.imag = z__2.imag + z__6.imag; // , expr subst
            v[2].real = z__1.real;
            v[2].imag = z__1.imag; // , expr subst
            i__1 = h_dim1 + 1;
            i__2 = h_dim1 * 3 + 3;
            z__5.real = h__[i__1].real + h__[i__2].real;
            z__5.imag = h__[i__1].imag + h__[i__2].imag; // , expr subst
            z__4.real = z__5.real - s1->real;
            z__4.imag = z__5.imag - s1->imag; // , expr subst
            z__3.real = z__4.real - s2->real;
            z__3.imag = z__4.imag - s2->imag; // , expr subst
            z__2.real = h31s.real * z__3.real - h31s.imag * z__3.imag;
            z__2.imag = h31s.real * z__3.imag + h31s.imag * z__3.real; // , expr subst
            i__3 = (h_dim1 << 1) + 3;
            z__6.real = h21s.real * h__[i__3].real - h21s.imag * h__[i__3].imag;
            z__6.imag = h21s.real * h__[i__3].imag + h21s.imag * h__[i__3].real; // , expr subst
            z__1.real = z__2.real + z__6.real;
            z__1.imag = z__2.imag + z__6.imag; // , expr subst
            v[3].real = z__1.real;
            v[3].imag = z__1.imag; // , expr subst
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* zlaqr1_ */
