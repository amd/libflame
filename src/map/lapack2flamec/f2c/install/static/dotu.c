#include "FLA_f2c.h"

#ifndef FLA_ENABLE_F2C_DOTC

extern
complex cdotu_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy);
VOID cdotu_f2c_(complex *r, integer *n, complex *cx, integer *incx, complex *cy, integer *incy)
{
    aocl_blas_cdotu(r, n, cx, incx, cy, incy);
}

void aocl_lapack_zdotu_f2c(dcomplex *r, aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy,
                aocl_int64_t *incy)
{
    aocl_blas_zdotu(r, n, cx, incx, cy, incy);
}

#endif
