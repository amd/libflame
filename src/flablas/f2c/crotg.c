/* crotg.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int crotg_(scomplex *ca, scomplex *cb, real *c__, scomplex *s)
{
    /* System generated locals */
    real r__1, r__2;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    double c_f2c_abs(scomplex *), sqrt(doublereal);
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    real norm;
    scomplex alpha;
    real scale;
    if (c_f2c_abs(ca) != 0.f)
    {
        goto L10;
    }
    *c__ = 0.f;
    s->real = 1.f, s->imag = 0.f;
    ca->real = cb->real, ca->imag = cb->imag;
    goto L20;
L10:
    scale = c_f2c_abs(ca) + c_f2c_abs(cb);
    q__1.real = ca->real / scale, q__1.imag = ca->imag / scale;
    /* Computing 2nd power */
    r__1 = c_f2c_abs(&q__1);
    q__2.real = cb->real / scale, q__2.imag = cb->imag / scale;
    /* Computing 2nd power */
    r__2 = c_f2c_abs(&q__2);
    norm = scale * sqrt(r__1 * r__1 + r__2 * r__2);
    r__1 = c_f2c_abs(ca);
    q__1.real = ca->real / r__1, q__1.imag = ca->imag / r__1;
    alpha.real = q__1.real, alpha.imag = q__1.imag;
    *c__ = c_f2c_abs(ca) / norm;
    r_cnjg(&q__3, cb);
    q__2.real = alpha.real * q__3.real - alpha.imag * q__3.imag, q__2.imag = alpha.real * q__3.imag + alpha.imag * q__3.real;
    q__1.real = q__2.real / norm, q__1.imag = q__2.imag / norm;
    s->real = q__1.real, s->imag = q__1.imag;
    q__1.real = norm * alpha.real, q__1.imag = norm * alpha.imag;
    ca->real = q__1.real, ca->imag = q__1.imag;
L20:
    return 0;
}
/* crotg_ */

