#include "FLA_f2c.h"
 void z_div(dcomplex *cp, dcomplex *ap, dcomplex *bp) {
 dcomplex a = *ap;
 dcomplex b = *bp;
 double temp;
 temp = b.r * b.r + b.i * b.i;
 cp->r = ( a.r * b.r + a.i * b.i ) / temp;
 cp->i = ( a.i * b.r - a.r * b.i ) / temp;
 }
 
