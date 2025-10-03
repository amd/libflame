#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

double d_imag(dcomplex *z)
{
    return (z->imag);
}
double r_imag(scomplex *z)
{
    return (z->imag);
}

#ifdef __cplusplus
}
#endif
