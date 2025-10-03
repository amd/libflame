#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif
/* Integer */
shortint h_abs(shortint *x)
{
    return (*x >= 0 ? (*x) : (-*x));
}
integer i_abs(integer *x)
{
    return (*x >= 0 ? (*x) : (-*x));
}

/* Double */
double r_abs(real *x)
{
    return (*x >= 0 ? (*x) : (-*x));
}
double d_abs(doublereal *x)
{
    return (*x >= 0 ? (*x) : (-*x));
}

#ifdef _WIN32
/* Complex */
double c_abs(scomplex *z)
{
    _Fcomplex z_ = {z->real, z->imag};
    return (cabsf(z_));
}
double z_abs(dcomplex *z)
{
    _Dcomplex z_ = {z->real, z->imag};
    return (cabs(z_));
}
#else
/* Complex */
double c_abs(scomplex *z)
{
    return (cabs(z->real + I * z->imag));
}
double z_abs(dcomplex *z)
{
    return (cabs(z->real + I * z->imag));
}
#endif

#ifdef __cplusplus
}
#endif
