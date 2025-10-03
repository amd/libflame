/* ../netlib/claesy.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static aocl_int64_t c__2 = 2;
/* > \brief \b CLAESY computes the eigenvalues and eigenvectors of a 2-by-2 scomplex symmetric
 * matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAESY + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claesy.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claesy.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claesy.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 ) */
/* .. Scalar Arguments .. */
/* COMPLEX A, B, C, CS1, EVSCAL, RT1, RT2, SN1 */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix */
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
/* > A is COMPLEX */
/* > The ( 1, 1 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX */
/* > The ( 1, 2 ) element of input matrix. The ( 2, 1 ) element */
/* > is also given by B, since the 2-by-2 matrix is symmetric. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is COMPLEX */
/* > The ( 2, 2 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1 */
/* > \verbatim */
/* > RT1 is COMPLEX */
/* > The eigenvalue of larger modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] RT2 */
/* > \verbatim */
/* > RT2 is COMPLEX */
/* > The eigenvalue of smaller modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] EVSCAL */
/* > \verbatim */
/* > EVSCAL is COMPLEX */
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
/* > CS1 is COMPLEX */
/* > \endverbatim */
/* > */
/* > \param[out] SN1 */
/* > \verbatim */
/* > SN1 is COMPLEX */
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
/* > \ingroup complexSYauxiliary */
/* ===================================================================== */
/* Subroutine */
void claesy_(scomplex *a, scomplex *b, scomplex *c__, scomplex *rt1, scomplex *rt2, scomplex *evscal,
             scomplex *cs1, scomplex *sn1)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
    /* System generated locals */
    real r__1, r__2;
    scomplex q__1, q__2, q__3, q__4, q__5, q__6, q__7;
    /* Builtin functions */
    double c_abs(scomplex *);
    void pow_ci(scomplex *, scomplex *, aocl_int64_t *), c_sqrt(scomplex *, scomplex *),
        c_div(scomplex *, scomplex *, scomplex *);
    /* Local variables */
    scomplex s, t;
    real z__;
    scomplex tmp;
    real babs, tabs, evnorm;
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
    if(c_abs(b) == 0.f)
    {
        rt1->real = a->real, rt1->imag = a->imag;
        rt2->real = c__->real, rt2->imag = c__->imag;
        if(c_abs(rt1) < c_abs(rt2))
        {
            tmp.real = rt1->real;
            tmp.imag = rt1->imag; // , expr subst
            rt1->real = rt2->real, rt1->imag = rt2->imag;
            rt2->real = tmp.real, rt2->imag = tmp.imag;
            cs1->real = 0.f, cs1->imag = 0.f;
            sn1->real = 1.f, sn1->imag = 0.f;
        }
        else
        {
            cs1->real = 1.f, cs1->imag = 0.f;
            sn1->real = 0.f, sn1->imag = 0.f;
        }
    }
    else
    {
        /* Compute the eigenvalues and eigenvectors. */
        /* The characteristic equation is */
        /* lambda **2 - (A+C) lambda + (A*C - B*B) */
        /* and we solve it using the quadratic formula. */
        q__2.real = a->real + c__->real;
        q__2.imag = a->imag + c__->imag; // , expr subst
        q__1.real = q__2.real * .5f;
        q__1.imag = q__2.imag * .5f; // , expr subst
        s.real = q__1.real;
        s.imag = q__1.imag; // , expr subst
        q__2.real = a->real - c__->real;
        q__2.imag = a->imag - c__->imag; // , expr subst
        q__1.real = q__2.real * .5f;
        q__1.imag = q__2.imag * .5f; // , expr subst
        t.real = q__1.real;
        t.imag = q__1.imag; // , expr subst
        /* Take the square root carefully to avoid over/under flow. */
        babs = c_abs(b);
        tabs = c_abs(&t);
        z__ = fla_max(babs, tabs);
        if(z__ > 0.f)
        {
            q__5.real = t.real / z__;
            q__5.imag = t.imag / z__; // , expr subst
            pow_ci(&q__4, &q__5, &c__2);
            q__7.real = b->real / z__;
            q__7.imag = b->imag / z__; // , expr subst
            pow_ci(&q__6, &q__7, &c__2);
            q__3.real = q__4.real + q__6.real;
            q__3.imag = q__4.imag + q__6.imag; // , expr subst
            c_sqrt(&q__2, &q__3);
            q__1.real = z__ * q__2.real;
            q__1.imag = z__ * q__2.imag; // , expr subst
            t.real = q__1.real;
            t.imag = q__1.imag; // , expr subst
        }
        /* Compute the two eigenvalues. RT1 and RT2 are exchanged */
        /* if necessary so that RT1 will have the greater magnitude. */
        q__1.real = s.real + t.real;
        q__1.imag = s.imag + t.imag; // , expr subst
        rt1->real = q__1.real, rt1->imag = q__1.imag;
        q__1.real = s.real - t.real;
        q__1.imag = s.imag - t.imag; // , expr subst
        rt2->real = q__1.real, rt2->imag = q__1.imag;
        if(c_abs(rt1) < c_abs(rt2))
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
        q__2.real = rt1->real - a->real;
        q__2.imag = rt1->imag - a->imag; // , expr subst
        c_div(&q__1, &q__2, b);
        sn1->real = q__1.real, sn1->imag = q__1.imag;
        tabs = c_abs(sn1);
        if(tabs > 1.f)
        {
            /* Computing 2nd power */
            r__2 = 1.f / tabs;
            r__1 = r__2 * r__2;
            q__5.real = sn1->real / tabs;
            q__5.imag = sn1->imag / tabs; // , expr subst
            pow_ci(&q__4, &q__5, &c__2);
            q__3.real = r__1 + q__4.real;
            q__3.imag = q__4.imag; // , expr subst
            c_sqrt(&q__2, &q__3);
            q__1.real = tabs * q__2.real;
            q__1.imag = tabs * q__2.imag; // , expr subst
            t.real = q__1.real;
            t.imag = q__1.imag; // , expr subst
        }
        else
        {
            q__3.real = sn1->real * sn1->real - sn1->imag * sn1->imag;
            q__3.imag = sn1->real * sn1->imag + sn1->imag * sn1->real; // , expr subst
            q__2.real = q__3.real + 1.f;
            q__2.imag = q__3.imag + 0.f; // , expr subst
            c_sqrt(&q__1, &q__2);
            t.real = q__1.real;
            t.imag = q__1.imag; // , expr subst
        }
        evnorm = c_abs(&t);
        if(evnorm >= .1f)
        {
            c_div(&q__1, &c_b1, &t);
            evscal->real = q__1.real, evscal->imag = q__1.imag;
            cs1->real = evscal->real, cs1->imag = evscal->imag;
            q__1.real = sn1->real * evscal->real - sn1->imag * evscal->imag;
            q__1.imag = sn1->real * evscal->imag + sn1->imag * evscal->real; // , expr subst
            sn1->real = q__1.real, sn1->imag = q__1.imag;
        }
        else
        {
            evscal->real = 0.f, evscal->imag = 0.f;
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLAESY */
}
/* claesy_ */
