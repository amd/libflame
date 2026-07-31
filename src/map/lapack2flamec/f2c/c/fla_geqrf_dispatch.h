/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_geqrf_dispatch.h
 *  @brief GEQRF multithreaded large path dispatch macros.
 *  */
#ifndef FLA_GEQRF_DISPATCH_H
#define FLA_GEQRF_DISPATCH_H

#include "FLAME.h"

#define FLA_SGEQRF_MT_EXTRA_OK(M, N, NUM_THREADS) \
    ((aocl_int64_t)(M) * (N) >= FLA_SGEQRF_MT_MIN_ELEMS)

#define FLA_CGEQRF_MT_EXTRA_OK(M, N, NUM_THREADS)                          \
    (((aocl_int64_t)(M) * (N) <= FLA_CGEQRF_MT_HUGE_THRESH                 \
      || (NUM_THREADS) >= FLA_GEQRF_MT_HUGE_MIN_THREADS)                   \
     && ((aocl_int64_t)(N) <= FLA_GEQRF_MT_WIDE_ASPECT * (aocl_int64_t)(M) \
         || (aocl_int64_t)(M) * (N) <= FLA_GEQRF_MT_WIDE_MIN_ELEMS)        \
     && !((aocl_int64_t)(M) * (N) > FLA_CGEQRF_MT_HUGE_THRESH && (N) > (M)))

#define FLA_ZGEQRF_MT_EXTRA_OK(M, N, NUM_THREADS)                          \
    (((aocl_int64_t)(M) * (N) <= FLA_ZGEQRF_MT_HUGE_THRESH                 \
      || (NUM_THREADS) >= FLA_GEQRF_MT_HUGE_MIN_THREADS)                   \
     && ((aocl_int64_t)(N) <= FLA_GEQRF_MT_WIDE_ASPECT * (aocl_int64_t)(M) \
         || (aocl_int64_t)(M) * (N) <= FLA_GEQRF_MT_WIDE_MIN_ELEMS))

#define FLA_GEQRF_MT_LARGE_DISPATCH(PRE, M, N, NUM_THREADS, LWORK, LWKOPT_MT)                      \
    ((NUM_THREADS) > 1 && (M) >= FLA_GEQRF_MT_LARGE_M_THRESH && (N) >= FLA_GEQRF_MT_LARGE_N_THRESH \
     && (LWORK) >= (LWKOPT_MT) && FLA_##PRE##GEQRF_MT_EXTRA_OK(M, N, NUM_THREADS))

#endif /* FLA_GEQRF_DISPATCH_H */
