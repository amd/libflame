#include "FLA_f2c.h"

#ifndef FLA_ENABLE_F2C_DOTC

extern
complex cdotc_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy);

void aocl_lapack_cdotc_f2c(scomplex *r, aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx, scomplex *cy, aocl_int64_t *incy)
{
    aocl_blas_cdotc(r, n, cx, incx, cy, incy);
}

void aocl_lapack_zdotc_f2c(dcomplex *r, aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy,
                aocl_int64_t *incy)
{
    aocl_blas_zdotc(r, n, cx, incx, cy, incy);
}
#endif
