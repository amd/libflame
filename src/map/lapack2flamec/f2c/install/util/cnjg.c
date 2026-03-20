#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

void d_cnjg(dcomplex *r, dcomplex *z)
{
    doublereal zi = z->imag;
    r->real = z->real;
    r->imag = -zi;
}
void r_cnjg(scomplex *r, scomplex *z)
{
    real zi = z->imag;
    r->real = z->real;
    r->imag = -zi;
}

#ifdef __cplusplus
}
#endif
