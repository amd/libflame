#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

void d_cnjg(dcomplex *r, dcomplex *z)
{
    doublereal zi = z->i;
    r->r = z->r;
    r->i = -zi;
}
void r_cnjg(scomplex *r, scomplex *z)
{
    real zi = z->i;
    r->r = z->r;
    r->i = -zi;
}

#ifdef __cplusplus
}
#endif
