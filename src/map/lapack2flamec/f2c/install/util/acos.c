#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

double r_acos(real *x)
{
    return (acos(*x));
}
double d_acos(doublereal *x)
{
    return (acos(*x));
}
/*
void c_acos(scomplex *r, scomplex *z)
{
  double _Complex ret_val = cacos(*z);
  r->real = creal(ret_val);
  r->imag = cimag(ret_val);
}
void z_acos(dcomplex *r, dcomplex *z)
{
  double _Complex ret_val = cacos(*z);
  r->real = creal(ret_val);
  r->imag = cimag(ret_val);
}
*/
#ifdef __cplusplus
}
#endif
