#include "FLA_f2c.h"

void aocl_lapack_cdotu_f2c(scomplex *r, aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx, scomplex *cy, aocl_int64_t *incy)
{
    aocl_blas_cdotu(r, n, cx, incx, cy, incy);
}

void aocl_lapack_zdotu_f2c(dcomplex *r, aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy,
                aocl_int64_t *incy)
{
    aocl_blas_zdotu(r, n, cx, incx, cy, incy);
}
