#include "FLA_f2c.h"
 void r_cnjg(scomplex *dest, scomplex *src) {
 dest->r = src->r ;
 dest->i = -(src->i);
 }
 
