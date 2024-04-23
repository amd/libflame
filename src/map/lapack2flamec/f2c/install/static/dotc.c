#include "FLA_f2c.h"

#ifndef FLA_ENABLE_F2C_DOTC

extern complex cdotc_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy);

void aocl_lapack_cdotc_f2c(scomplex *r, aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx, scomplex *cy, aocl_int64_t *incy)
{
    *r = cdotc_(n, cx, incx, cy, incy);
}

extern doublecomplex zdotc_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy,
                            integer *incy);
VOID zdotc_f2c_(doublecomplex *r, integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy,
                integer *incy)
{
    *r = zdotc_(n, cx, incx, cy, incy);
}
#endif
