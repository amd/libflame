#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

double r_cos(real *x)
{
    return (cos(*x));
}
double d_cos(doublereal *x)
{
    return (cos(*x));
}

#ifdef _WIN32
void c_cos(scomplex *r, scomplex *z)
{
    _Fcomplex z_ = {z->real, z->imag};
    _Fcomplex ret_val = ccosf(z_);
    r->real = crealf(ret_val);
    r->imag = cimagf(ret_val);
}
void z_cos(dcomplex *r, dcomplex *z)
{
    _Dcomplex z_ = {z->real, z->imag};
    _Dcomplex ret_val = ccos(z_);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#else
void c_cos(scomplex *r, scomplex *z)
{
    double _Complex ret_val = ccos(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
void z_cos(dcomplex *r, dcomplex *z)
{
    double _Complex ret_val = ccos(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#endif

#ifdef __cplusplus
}
#endif
