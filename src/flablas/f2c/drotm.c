/* drotm.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int drotm_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *dparam)
{
    /* Initialized data */
    static const doublereal zero = 0.;
    static const doublereal two = 2.;
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    integer i__;
    doublereal dflag, w, z__;
    integer kx, ky, nsteps;
    doublereal dh11, dh12, dh22, dh21;
    /* APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */
    /* (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN */
    /* (DY**T) */
    /* DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE */
    /* LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY. */
    /* WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
    /* DFLAG=-1.D0 DFLAG=0.D0 DFLAG=1.D0 DFLAG=-2.D0 */
    /* (DH11 DH12) (1.D0 DH12) (DH11 1.D0) (1.D0 0.D0) */
    /* H=( ) ( ) ( ) ( ) */
    /* (DH21 DH22), (DH21 1.D0), (-1.D0 DH22), (0.D0 1.D0). */
    /* SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM. */
    /* Parameter adjustments */
    --dparam;
    --dy;
    --dx;
    /* Function Body */
    dflag = dparam[1];
    if (*n <= 0 || dflag + two == zero)
    {
        goto L140;
    }
    if (! (*incx == *incy && *incx > 0))
    {
        goto L70;
    }
    nsteps = *n * *incx;
    if (dflag < 0.)
    {
        goto L50;
    }
    else if (dflag == 0)
    {
        goto L10;
    }
    else
    {
        goto L30;
    }
L10:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        w = dx[i__];
        z__ = dy[i__];
        dx[i__] = w + z__ * dh12;
        dy[i__] = w * dh21 + z__;
        /* L20: */
    }
    goto L140;
L30:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = nsteps;
    i__1 = *incx;
    for (i__ = 1;
            i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
            i__ += i__1)
    {
        w = dx[i__];
        z__ = dy[i__];
        dx[i__] = w * dh11 + z__;
        dy[i__] = -w + dh22 * z__;
        /* L40: */
    }
    goto L140;
L50:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        w = dx[i__];
        z__ = dy[i__];
        dx[i__] = w * dh11 + z__ * dh12;
        dy[i__] = w * dh21 + z__ * dh22;
        /* L60: */
    }
    goto L140;
L70:
    kx = 1;
    ky = 1;
    if (*incx < 0)
    {
        kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0)
    {
        ky = (1 - *n) * *incy + 1;
    }
    if (dflag < 0.)
    {
        goto L120;
    }
    else if (dflag == 0)
    {
        goto L80;
    }
    else
    {
        goto L100;
    }
L80:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__2 = *n;
    for (i__ = 1;
            i__ <= i__2;
            ++i__)
    {
        w = dx[kx];
        z__ = dy[ky];
        dx[kx] = w + z__ * dh12;
        dy[ky] = w * dh21 + z__;
        kx += *incx;
        ky += *incy;
        /* L90: */
    }
    goto L140;
L100:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1;
            i__ <= i__2;
            ++i__)
    {
        w = dx[kx];
        z__ = dy[ky];
        dx[kx] = w * dh11 + z__;
        dy[ky] = -w + dh22 * z__;
        kx += *incx;
        ky += *incy;
        /* L110: */
    }
    goto L140;
L120:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1;
            i__ <= i__2;
            ++i__)
    {
        w = dx[kx];
        z__ = dy[ky];
        dx[kx] = w * dh11 + z__ * dh12;
        dy[ky] = w * dh21 + z__ * dh22;
        kx += *incx;
        ky += *incy;
        /* L130: */
    }
L140:
    return 0;
}
/* drotm_ */

