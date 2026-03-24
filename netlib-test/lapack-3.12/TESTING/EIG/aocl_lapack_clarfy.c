/**
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

#include <stdint.h>

#ifdef INT_64BIT
typedef int64_t aocl_int_t;
#else
typedef int32_t aocl_int_t;
#endif

typedef struct scomplex_
{
    float real;
    float imag;
} scomplex;

void aocl_lapack_clarfy(char *uplo, int64_t *n, scomplex *v, int64_t *incv, scomplex *tau,
                        scomplex *c__, int64_t *ldc, scomplex *work)
{
    /* System generated locals */
    extern void clarfy_(char *, aocl_int_t *, scomplex *, aocl_int_t *, scomplex *, scomplex *,
                        aocl_int_t *, scomplex *);

    aocl_int_t n_lp = (aocl_int_t)(*n);
    aocl_int_t incv_lp = (aocl_int_t)(*incv);
    aocl_int_t ldc_lp = (aocl_int_t)(*ldc);

    clarfy_(uplo, &n_lp, v, &incv_lp, tau, c__, &ldc_lp, work);
}
