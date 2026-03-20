#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

double r_sqrt(real *x)
{
    return (sqrt(*x));
}
double d_sqrt(doublereal *x)
{
    return (sqrt(*x));
}
#ifdef _WIN32
void c_sqrt(scomplex *r, scomplex *z)
{
    _Fcomplex z_ = {z->real, z->imag};
    _Fcomplex ret_val = csqrtf(z_);
    r->real = crealf(ret_val);
    r->imag = cimagf(ret_val);
}
void z_sqrt(dcomplex *r, dcomplex *z)
{
    _Dcomplex z_ = {z->real, z->imag};
    _Dcomplex ret_val = csqrt(z_);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#else
void c_sqrt(scomplex *r, scomplex *z)
{
    double _Complex ret_val = csqrt(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
void z_sqrt(dcomplex *r, dcomplex *z)
{
    double _Complex ret_val = csqrt(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#endif

#ifdef __cplusplus
}
#endif
