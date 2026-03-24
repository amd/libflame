/**
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

#include <stdint.h>

#ifdef INT_64BIT
typedef int64_t aocl_int_t;
#else
typedef int32_t aocl_int_t;
#endif

typedef double doublereal;

void aocl_lapack_dlarfy(char *uplo, int64_t *n, doublereal *v, int64_t *incv, doublereal *tau,
                        doublereal *c__, int64_t *ldc, doublereal *work)
{
    /* System generated locals */
    extern void dlarfy_(char *, aocl_int_t *, doublereal *, aocl_int_t *, doublereal *,
                        doublereal *, aocl_int_t *, doublereal *);

    aocl_int_t n_lp = (aocl_int_t)(*n);
    aocl_int_t incv_lp = (aocl_int_t)(*incv);
    aocl_int_t ldc_lp = (aocl_int_t)(*ldc);

    dlarfy_(uplo, &n_lp, v, &incv_lp, tau, c__, &ldc_lp, work);
}
