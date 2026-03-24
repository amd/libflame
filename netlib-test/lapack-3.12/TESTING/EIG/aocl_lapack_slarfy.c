/**
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

#include <stdint.h>

#ifdef INT_64BIT
typedef int64_t aocl_int_t;
#else
typedef int32_t aocl_int_t;
#endif

typedef float real;

void aocl_lapack_slarfy(char *uplo, int64_t *n, real *v, int64_t *incv, real *tau, real *c__,
                        int64_t *ldc, real *work)
{
    /* System generated locals */
    extern void slarfy_(char *, aocl_int_t *, real *, aocl_int_t *, real *, real *, aocl_int_t *,
                        real *);

    aocl_int_t n_lp = (aocl_int_t)(*n);
    aocl_int_t incv_lp = (aocl_int_t)(*incv);
    aocl_int_t ldc_lp = (aocl_int_t)(*ldc);

    slarfy_(uplo, &n_lp, v, &incv_lp, tau, c__, &ldc_lp, work);
}
