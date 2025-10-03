#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

double r_log(real *x)
{
    return (log(*x));
}
double d_log(doublereal *x)
{
    return (log(*x));
}

#ifdef _WIN32
void c_log(scomplex *r, scomplex *z)
{
    _Fcomplex z_ = {z->real, z->imag};
    _Fcomplex ret_val = clogf(z_);
    r->real = crealf(ret_val);
    r->imag = cimagf(ret_val);
}
void z_log(dcomplex *r, dcomplex *z)
{
    _Dcomplex z_ = {z->real, z->imag};
    _Dcomplex ret_val = clog(z_);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#else
void c_log(scomplex *r, scomplex *z)
{
    double _Complex ret_val = clog(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
void z_log(dcomplex *r, dcomplex *z)
{
    double _Complex ret_val = clog(z->real + I * z->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
#endif

#ifdef __cplusplus
}
#endif
