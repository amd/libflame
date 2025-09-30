#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

double d_imag(dcomplex *z)
{
    return (z->i);
}
double r_imag(scomplex *z)
{
    return (z->i);
}

#ifdef __cplusplus
}
#endif
