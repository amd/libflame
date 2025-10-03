#include "FLA_f2c.h"
 void d_cnjg(dcomplex *dest, dcomplex *src) {
 dest->real = src->real ;
 dest->imag = -(src->imag);
 }
 
