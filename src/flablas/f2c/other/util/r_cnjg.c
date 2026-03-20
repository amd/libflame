#include "FLA_f2c.h"
 void r_cnjg(scomplex *dest, scomplex *src) {
 dest->real = src->real ;
 dest->imag = -(src->imag);
 }
 
