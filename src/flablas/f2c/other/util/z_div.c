#include "FLA_f2c.h"
 void z_div(dcomplex *cp, dcomplex *ap, dcomplex *bp) {
 dcomplex a = *ap;
 dcomplex b = *bp;
 double temp;
 temp = b.real * b.real + b.imag * b.imag;
 cp->real = ( a.real * b.real + a.imag * b.imag ) / temp;
 cp->imag = ( a.imag * b.real - a.real * b.imag ) / temp;
 }
 
