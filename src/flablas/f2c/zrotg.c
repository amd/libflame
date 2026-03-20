/* zrotg.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zrotg_(dcomplex *ca, dcomplex *cb, doublereal * c__, dcomplex *s)
{
    /* System generated locals */
    doublereal d__1, d__2;
    dcomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    double z_f2c_abs(dcomplex *);
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    double sqrt(doublereal);
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    doublereal norm;
    dcomplex alpha;
    doublereal scale;
    if (z_f2c_abs(ca) != 0.)
    {
        goto L10;
    }
    *c__ = 0.;
    s->real = 1., s->imag = 0.;
    ca->real = cb->real, ca->imag = cb->imag;
    goto L20;
L10:
    scale = z_f2c_abs(ca) + z_f2c_abs(cb);
    z__2.real = scale, z__2.imag = 0.;
    z_div(&z__1, ca, &z__2);
    /* Computing 2nd power */
    d__1 = z_f2c_abs(&z__1);
    z__4.real = scale, z__4.imag = 0.;
    z_div(&z__3, cb, &z__4);
    /* Computing 2nd power */
    d__2 = z_f2c_abs(&z__3);
    norm = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    d__1 = z_f2c_abs(ca);
    z__1.real = ca->real / d__1, z__1.imag = ca->imag / d__1;
    alpha.real = z__1.real, alpha.imag = z__1.imag;
    *c__ = z_f2c_abs(ca) / norm;
    d_cnjg(&z__3, cb);
    z__2.real = alpha.real * z__3.real - alpha.imag * z__3.imag, z__2.imag = alpha.real * z__3.imag + alpha.imag * z__3.real;
    z__1.real = z__2.real / norm, z__1.imag = z__2.imag / norm;
    s->real = z__1.real, s->imag = z__1.imag;
    z__1.real = norm * alpha.real, z__1.imag = norm * alpha.imag;
    ca->real = z__1.real, ca->imag = z__1.imag;
L20:
    return 0;
}
/* zrotg_ */

