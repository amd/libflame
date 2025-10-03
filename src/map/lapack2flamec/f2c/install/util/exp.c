#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

double r_exp(real *x)
{
    return (exp(*x));
}
double d_exp(doublereal *x)
{
    return (exp(*x));
}

#ifdef _WIN32
void c_exp(scomplex *r, scomplex *z)
{
    _Fcomplex z_ = {z->real, z->imag};
    _Fcomplex ret_val = cexpf(z_);
    r->real = crealf(ret_val);
    r->imag = cimagf(ret_val);
}
void z_exp(dcomplex *r, dcomplex *z)
{
    _Dcomplex z_ = {z->real, z->imag};
    _Dcomplex ret_val = cexp(z_);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#else
void c_exp(scomplex *r, scomplex *z)
{
    double _Complex ret_val = cexp(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
void z_exp(dcomplex *r, dcomplex *z)
{
    double _Complex ret_val = cexp(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#endif

#ifdef __cplusplus
}
#endif
