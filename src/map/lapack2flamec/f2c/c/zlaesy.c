/* ../netlib/zlaesy.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {1., 0.};
static aocl_int64_t c__2 = 2;
/* > \brief \b ZLAESY computes the eigenvalues and eigenvectors of a 2-by-2 scomplex symmetric
 * matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAESY + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaesy.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaesy.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaesy.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 ) */
/* .. Scalar Arguments .. */
/* COMPLEX*16 A, B, C, CS1, EVSCAL, RT1, RT2, SN1 */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix */
/* > ( ( A, B );
( B, C ) ) */
/* > provided the norm of the matrix of eigenvectors is larger than */
/* > some threshold value. */
/* > */
/* > RT1 is the eigenvalue of larger absolute value, and RT2 of */
/* > smaller absolute value. If the eigenvectors are computed, then */
/* > on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence */
/* > */
/* > [ CS1 SN1 ] . [ A B ] . [ CS1 -SN1 ] = [ RT1 0 ] */
/* > [ -SN1 CS1 ] [ B C ] [ SN1 CS1 ] [ 0 RT2 ] */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 */
/* > The ( 1, 1 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX*16 */
/* > The ( 1, 2 ) element of input matrix. The ( 2, 1 ) element */
/* > is also given by B, since the 2-by-2 matrix is symmetric. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is COMPLEX*16 */
/* > The ( 2, 2 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1 */
/* > \verbatim */
/* > RT1 is COMPLEX*16 */
/* > The eigenvalue of larger modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] RT2 */
/* > \verbatim */
/* > RT2 is COMPLEX*16 */
/* > The eigenvalue of smaller modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] EVSCAL */
/* > \verbatim */
/* > EVSCAL is COMPLEX*16 */
/* > The scomplex value by which the eigenvector matrix was scaled */
/* > to make it orthonormal. If EVSCAL is zero, the eigenvectors */
/* > were not computed. This means one of two things: the 2-by-2 */
/* > matrix could not be diagonalized, or the norm of the matrix */
/* > of eigenvectors before scaling was larger than the threshold */
/* > value THRESH (set below). */
/* > \endverbatim */
/* > */
/* > \param[out] CS1 */
/* > \verbatim */
/* > CS1 is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[out] SN1 */
/* > \verbatim */
/* > SN1 is COMPLEX*16 */
/* > If EVSCAL .NE. 0, ( CS1, SN1 ) is the unit right eigenvector */
/* > for RT1. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16SYauxiliary */
/* ===================================================================== */
/* Subroutine */
void zlaesy_(dcomplex *a, dcomplex *b, dcomplex *c__, dcomplex *rt1,
             dcomplex *rt2, dcomplex *evscal, dcomplex *cs1, dcomplex *sn1)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    /* System generated locals */
    doublereal d__1, d__2;
    dcomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;
    /* Builtin functions */
    double z_abs(dcomplex *);
    void pow_zi(dcomplex *, dcomplex *, aocl_int64_t *),
        z_sqrt(dcomplex *, dcomplex *),
        z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    dcomplex s, t;
    doublereal z__;
    dcomplex tmp;
    doublereal babs, tabs, evnorm;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Special case: The matrix is actually diagonal. */
    /* To avoid divide by zero later, we treat this case separately. */
    if(z_abs(b) == 0.)
    {
        rt1->real = a->real, rt1->imag = a->imag;
        rt2->real = c__->real, rt2->imag = c__->imag;
        if(z_abs(rt1) < z_abs(rt2))
        {
            tmp.real = rt1->real;
            tmp.imag = rt1->imag; // , expr subst
            rt1->real = rt2->real, rt1->imag = rt2->imag;
            rt2->real = tmp.real, rt2->imag = tmp.imag;
            cs1->real = 0., cs1->imag = 0.;
            sn1->real = 1., sn1->imag = 0.;
        }
        else
        {
            cs1->real = 1., cs1->imag = 0.;
            sn1->real = 0., sn1->imag = 0.;
        }
    }
    else
    {
        /* Compute the eigenvalues and eigenvectors. */
        /* The characteristic equation is */
        /* lambda **2 - (A+C) lambda + (A*C - B*B) */
        /* and we solve it using the quadratic formula. */
        z__2.real = a->real + c__->real;
        z__2.imag = a->imag + c__->imag; // , expr subst
        z__1.real = z__2.real * .5;
        z__1.imag = z__2.imag * .5; // , expr subst
        s.real = z__1.real;
        s.imag = z__1.imag; // , expr subst
        z__2.real = a->real - c__->real;
        z__2.imag = a->imag - c__->imag; // , expr subst
        z__1.real = z__2.real * .5;
        z__1.imag = z__2.imag * .5; // , expr subst
        t.real = z__1.real;
        t.imag = z__1.imag; // , expr subst
        /* Take the square root carefully to avoid over/under flow. */
        babs = z_abs(b);
        tabs = z_abs(&t);
        z__ = fla_max(babs, tabs);
        if(z__ > 0.)
        {
            z__5.real = t.real / z__;
            z__5.imag = t.imag / z__; // , expr subst
            pow_zi(&z__4, &z__5, &c__2);
            z__7.real = b->real / z__;
            z__7.imag = b->imag / z__; // , expr subst
            pow_zi(&z__6, &z__7, &c__2);
            z__3.real = z__4.real + z__6.real;
            z__3.imag = z__4.imag + z__6.imag; // , expr subst
            z_sqrt(&z__2, &z__3);
            z__1.real = z__ * z__2.real;
            z__1.imag = z__ * z__2.imag; // , expr subst
            t.real = z__1.real;
            t.imag = z__1.imag; // , expr subst
        }
        /* Compute the two eigenvalues. RT1 and RT2 are exchanged */
        /* if necessary so that RT1 will have the greater magnitude. */
        z__1.real = s.real + t.real;
        z__1.imag = s.imag + t.imag; // , expr subst
        rt1->real = z__1.real, rt1->imag = z__1.imag;
        z__1.real = s.real - t.real;
        z__1.imag = s.imag - t.imag; // , expr subst
        rt2->real = z__1.real, rt2->imag = z__1.imag;
        if(z_abs(rt1) < z_abs(rt2))
        {
            tmp.real = rt1->real;
            tmp.imag = rt1->imag; // , expr subst
            rt1->real = rt2->real, rt1->imag = rt2->imag;
            rt2->real = tmp.real, rt2->imag = tmp.imag;
        }
        /* Choose CS1 = 1 and SN1 to satisfy the first equation, then */
        /* scale the components of this eigenvector so that the matrix */
        /* of eigenvectors X satisfies X * X**T = I . (No scaling is */
        /* done if the norm of the eigenvalue matrix is less than THRESH.) */
        z__2.real = rt1->real - a->real;
        z__2.imag = rt1->imag - a->imag; // , expr subst
        z_div(&z__1, &z__2, b);
        sn1->real = z__1.real, sn1->imag = z__1.imag;
        tabs = z_abs(sn1);
        if(tabs > 1.)
        {
            /* Computing 2nd power */
            d__2 = 1. / tabs;
            d__1 = d__2 * d__2;
            z__5.real = sn1->real / tabs;
            z__5.imag = sn1->imag / tabs; // , expr subst
            pow_zi(&z__4, &z__5, &c__2);
            z__3.real = d__1 + z__4.real;
            z__3.imag = z__4.imag; // , expr subst
            z_sqrt(&z__2, &z__3);
            z__1.real = tabs * z__2.real;
            z__1.imag = tabs * z__2.imag; // , expr subst
            t.real = z__1.real;
            t.imag = z__1.imag; // , expr subst
        }
        else
        {
            z__3.real = sn1->real * sn1->real - sn1->imag * sn1->imag;
            z__3.imag = sn1->real * sn1->imag + sn1->imag * sn1->real; // , expr subst
            z__2.real = z__3.real + 1.;
            z__2.imag = z__3.imag + 0.; // , expr subst
            z_sqrt(&z__1, &z__2);
            t.real = z__1.real;
            t.imag = z__1.imag; // , expr subst
        }
        evnorm = z_abs(&t);
        if(evnorm >= .1)
        {
            z_div(&z__1, &c_b1, &t);
            evscal->real = z__1.real, evscal->imag = z__1.imag;
            cs1->real = evscal->real, cs1->imag = evscal->imag;
            z__1.real = sn1->real * evscal->real - sn1->imag * evscal->imag;
            z__1.imag = sn1->real * evscal->imag + sn1->imag * evscal->real; // , expr subst
            sn1->real = z__1.real, sn1->imag = z__1.imag;
        }
        else
        {
            evscal->real = 0., evscal->imag = 0.;
        }
    }
    AOCL_DTL_TRACE_EXIT_INDENT
    return;
    /* End of ZLAESY */
}
/* zlaesy_ */
