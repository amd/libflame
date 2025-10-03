#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Integer */
shortint pow_hh(shortint *ap, shortint *bp)
{
    return (shortint)(pow(*ap, *bp));
}
integer pow_ii(integer *ap, integer *bp)
{
    return (integer)(pow(*ap, *bp));
}
#ifdef INTEGER_STAR_8
longint pow_qq(longint *ap, longint *bp)
{
    return (longint)(pow(*ap, *bp));
}
#endif

/* Double */
double pow_ri(real *ap, integer *bp)
{
    return (pow(*ap, *bp));
}
double pow_dd(doublereal *ap, doublereal *bp)
{
    return (pow(*ap, *bp));
}
double pow_di(doublereal *ap, integer *bp)
{
    return (pow(*ap, *bp));
}

#ifdef _WIN32
/* Complex */
void pow_ci(scomplex *p, scomplex *a, integer *b)
{
    _Fcomplex a_ = {a->real, a->imag};
    _Fcomplex b_ = {*b, 0};
    _Fcomplex ret_val = cpowf(a_, b_);
    p->real = crealf(ret_val);
    p->imag = cimagf(ret_val);
}
void pow_zz(dcomplex *r, dcomplex *a, dcomplex *b)
{
    _Dcomplex a_ = {a->real, a->imag};
    _Dcomplex b_ = {a->real, a->imag};
    _Dcomplex ret_val = cpow(a_, b_);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
void pow_zi(dcomplex *p, dcomplex *a, integer *b)
{
    _Dcomplex a_ = {a->real, a->imag};
    _Dcomplex b_ = {*b, 0};
    _Dcomplex ret_val = cpow(a_, b_);
    p->real = creal(ret_val);
    p->imag = cimag(ret_val);
}
#else
/* Complex */
void pow_ci(scomplex *p, scomplex *a, integer *b)
{
    double _Complex ret_val = cpow(a->real + I * a->imag, *b);
    p->real = creal(ret_val);
    p->imag = cimag(ret_val);
}
void pow_zz(dcomplex *r, dcomplex *a, dcomplex *b)
{
    double _Complex ret_val = cpow(a->real + I * a->imag, b->real + I * b->imag);
    r->real = creal(ret_val);
    r->imag = cimag(ret_val);
}
void pow_zi(dcomplex *p, dcomplex *a, integer *b)
{
    double _Complex ret_val = cpow(a->real + I * a->imag, *b);
    p->real = creal(ret_val);
    p->imag = cimag(ret_val);
}
#endif

#ifdef __cplusplus
}
#endif
