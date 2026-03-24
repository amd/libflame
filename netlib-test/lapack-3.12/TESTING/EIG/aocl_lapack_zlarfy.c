/**
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

#include <stdint.h>

#ifdef INT_64BIT
typedef int64_t aocl_int_t;
#else
typedef int32_t aocl_int_t;
#endif

typedef struct dcomplex_
{
    double real;
    double imag;
} dcomplex;

void aocl_lapack_zlarfy(char *uplo, int64_t *n, dcomplex *v, int64_t *incv, dcomplex *tau,
                        dcomplex *c__, int64_t *ldc, dcomplex *work)
{
    /* System generated locals */
    extern void zlarfy_(char *, aocl_int_t *, dcomplex *, aocl_int_t *, dcomplex *, dcomplex *,
                        aocl_int_t *, dcomplex *);

    aocl_int_t n_lp = (aocl_int_t)(*n);
    aocl_int_t incv_lp = (aocl_int_t)(*incv);
    aocl_int_t ldc_lp = (aocl_int_t)(*ldc);

    zlarfy_(uplo, &n_lp, v, &incv_lp, tau, c__, &ldc_lp, work);
}
