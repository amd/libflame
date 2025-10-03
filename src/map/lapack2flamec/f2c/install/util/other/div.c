#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
void c_div(scomplex *c, scomplex *a, scomplex *b)
{
    c->real = ((a->real * b->real) + (a->imag * b->imag)) / ((b->real * b->real) + (b->imag * b->imag));
    c->imag = ((a->imag * b->real) - (a->real * b->imag)) / ((b->real * b->real) + (b->imag * b->imag));
}
void z_div(dcomplex *c, dcomplex *a, dcomplex *b)
{
    c->real = ((a->real * b->real) + (a->imag * b->imag)) / ((b->real * b->real) + (b->imag * b->imag));
    c->imag = ((a->imag * b->real) - (a->real * b->imag)) / ((b->real * b->real) + (b->imag * b->imag));
}
#else
void c_div(scomplex *c, scomplex *a, scomplex *b)
{
    double _Complex ret_val = (a->real + I * a->imag) / (b->real + I * b->imag);
    c->real = creal(ret_val);
    c->imag = cimag(ret_val);
}
void z_div(dcomplex *c, dcomplex *a, dcomplex *b)
{
    double _Complex ret_val = (a->real + I * a->imag) / (b->real + I * b->imag);
    c->real = creal(ret_val);
    c->imag = cimag(ret_val);
}
#endif

#ifdef __cplusplus
}
#endif
