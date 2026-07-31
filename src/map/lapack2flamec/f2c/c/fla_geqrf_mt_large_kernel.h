/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*
 * Multithreaded large GEQRF. Set GEQRF_MT_PRECISION to D, S, C, or Z before
 * including this header.
 */

#if FLA_ENABLE_AMD_OPT && FLA_OPENMP_MULTITHREADING

#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

#ifndef FLA_GEQRF_MT_LARGE_KERNEL_BINDINGS
#define FLA_GEQRF_MT_LARGE_KERNEL_BINDINGS

#define GEQRF_MT_CONCAT_IMPL(a, b) a##b
#define GEQRF_MT_CONCAT(a, b) GEQRF_MT_CONCAT_IMPL(a, b)

/* Double (D) */
#define GEQRF_MT_TYPE_D doublereal
#define GEQRF_MT_GEQR2_D dgeqr2_fla
#define GEQRF_MT_LARFB_D aocl_lapack_dlarfb
#define GEQRF_MT_LARFT_D aocl_lapack_dlarft
#define GEQRF_MT_LARGE_D dgeqrf_mt_large
#define GEQRF_MT_THREAD_D dgeqrf_large_mt_thread_fn
#define GEQRF_MT_NUM_THREADS_D dgeqrf_mt_large_num_threads
#define GEQRF_MT_LWORK_D dgeqrf_mt_large_lwork
#define GEQRF_MT_PANEL_SIZE_D FLA_GEQRF_MT_LARGE_PANEL_SIZE
#define GEQRF_MT_THRESHOLD_8_THREADS_D FLA_GEQRF_MT_THRESHOLD_8_THREADS_REAL
#define GEQRF_MT_TRANS_D "Transpose"

/* Single (S) */
#define GEQRF_MT_TYPE_S real
#define GEQRF_MT_GEQR2_S sgeqr2_fla
#define GEQRF_MT_LARFB_S aocl_lapack_slarfb
#define GEQRF_MT_LARFT_S aocl_lapack_slarft
#define GEQRF_MT_LARGE_S sgeqrf_mt_large
#define GEQRF_MT_THREAD_S sgeqrf_large_mt_thread_fn
#define GEQRF_MT_NUM_THREADS_S sgeqrf_mt_large_num_threads
#define GEQRF_MT_LWORK_S sgeqrf_mt_large_lwork
#define GEQRF_MT_PANEL_SIZE_S FLA_GEQRF_MT_LARGE_PANEL_SIZE
#define GEQRF_MT_THRESHOLD_8_THREADS_S FLA_GEQRF_MT_THRESHOLD_8_THREADS_REAL
#define GEQRF_MT_TRANS_S "Transpose"

/* Complex single (C) */
#define GEQRF_MT_TYPE_C scomplex
#define GEQRF_MT_GEQR2_C aocl_lapack_cgeqr2
#define GEQRF_MT_LARFB_C aocl_lapack_clarfb
#define GEQRF_MT_LARFT_C aocl_lapack_clarft
#define GEQRF_MT_LARGE_C cgeqrf_mt_large
#define GEQRF_MT_THREAD_C cgeqrf_large_mt_thread_fn
#define GEQRF_MT_NUM_THREADS_C cgeqrf_mt_large_num_threads
#define GEQRF_MT_LWORK_C cgeqrf_mt_large_lwork
#define GEQRF_MT_PANEL_SIZE_C FLA_GEQRF_MT_LARGE_PANEL_SIZE
#define GEQRF_MT_THRESHOLD_8_THREADS_C FLA_CGEQRF_MT_THRESHOLD_8_THREADS
#define GEQRF_MT_TRANS_C "Conjugate transpose"

/* Complex double (Z) */
#define GEQRF_MT_TYPE_Z dcomplex
#define GEQRF_MT_GEQR2_Z aocl_lapack_zgeqr2
#define GEQRF_MT_LARFB_Z aocl_lapack_zlarfb
#define GEQRF_MT_LARFT_Z aocl_lapack_zlarft
#define GEQRF_MT_LARGE_Z zgeqrf_mt_large
#define GEQRF_MT_THREAD_Z zgeqrf_large_mt_thread_fn
#define GEQRF_MT_NUM_THREADS_Z zgeqrf_mt_large_num_threads
#define GEQRF_MT_LWORK_Z zgeqrf_mt_large_lwork
#define GEQRF_MT_PANEL_SIZE_Z FLA_GEQRF_MT_LARGE_PANEL_SIZE
#define GEQRF_MT_THRESHOLD_8_THREADS_Z FLA_ZGEQRF_MT_THRESHOLD_8_THREADS
#define GEQRF_MT_TRANS_Z "Conjugate transpose"

static inline size_t containers_needed_for_size(aocl_int64_t nb, size_t item_size,
                                                size_t container_size)
{
    size_t size_needed = item_size * nb;
    return (size_needed - 1) / container_size + 1;
}

#endif /* FLA_GEQRF_MT_LARGE_KERNEL_BINDINGS */

#ifndef GEQRF_MT_PRECISION
#error "fla_geqrf_mt_large_kernel.h: define GEQRF_MT_PRECISION to D, S, C, or Z before including"
#endif

#define GEQRF_MT_TYPE GEQRF_MT_CONCAT(GEQRF_MT_TYPE_, GEQRF_MT_PRECISION)
#define GEQRF_MT_GEQR2_FN GEQRF_MT_CONCAT(GEQRF_MT_GEQR2_, GEQRF_MT_PRECISION)
#define GEQRF_MT_LARFB_FN GEQRF_MT_CONCAT(GEQRF_MT_LARFB_, GEQRF_MT_PRECISION)
#define GEQRF_MT_LARFT_FN GEQRF_MT_CONCAT(GEQRF_MT_LARFT_, GEQRF_MT_PRECISION)
#define GEQRF_MT_LARGE_FN GEQRF_MT_CONCAT(GEQRF_MT_LARGE_, GEQRF_MT_PRECISION)
#define GEQRF_MT_THREAD_FN GEQRF_MT_CONCAT(GEQRF_MT_THREAD_, GEQRF_MT_PRECISION)
#define GEQRF_MT_NUM_THREADS_FN GEQRF_MT_CONCAT(GEQRF_MT_NUM_THREADS_, GEQRF_MT_PRECISION)
#define GEQRF_MT_LWORK_FN GEQRF_MT_CONCAT(GEQRF_MT_LWORK_, GEQRF_MT_PRECISION)
#define GEQRF_MT_PANEL_SIZE GEQRF_MT_CONCAT(GEQRF_MT_PANEL_SIZE_, GEQRF_MT_PRECISION)
#define GEQRF_MT_THRESHOLD_8_THREADS \
    GEQRF_MT_CONCAT(GEQRF_MT_THRESHOLD_8_THREADS_, GEQRF_MT_PRECISION)
#define GEQRF_MT_TRANS GEQRF_MT_CONCAT(GEQRF_MT_TRANS_, GEQRF_MT_PRECISION)

/* Returns the optimial number of threads to use for the mt_large path */
static integer GEQRF_MT_NUM_THREADS_FN(aocl_int64_t gm, aocl_int64_t gn)
{
    aocl_int64_t num_elems = gm * gn;
    aocl_int64_t opt_nthreads;
    if(num_elems <= GEQRF_MT_THRESHOLD_8_THREADS)
    {
        opt_nthreads = 8;
    }
    else
    {
        opt_nthreads = 64;
    }
    return fla_min(fla_thread_get_num_threads(), opt_nthreads);
}

/* Returns the optimal lwork for multithreaded geqrf */
static integer GEQRF_MT_LWORK_FN(aocl_int64_t gm, aocl_int64_t gn, aocl_int64_t num_threads)
{

    aocl_int64_t nb = GEQRF_MT_PANEL_SIZE;

    /* Number of panels */
    aocl_int64_t nt = ((gn - 1) / nb) + 1;

    /* Storage to keep triangular factors for each panel */
    aocl_int64_t triangular_factor_req = (nb * nb) * nt;

    aocl_int64_t per_thread_work_req = (nb * nb) * num_threads;

    /* Dependency list to keep track of which panels are ready */
    aocl_int64_t depend_list_req
        = containers_needed_for_size(nt, sizeof(uint8_t), sizeof(GEQRF_MT_TYPE));

    aocl_int64_t lwork = triangular_factor_req + per_thread_work_req + depend_list_req;

    return lwork;
}

static void GEQRF_MT_THREAD_FN(aocl_int64_t m, aocl_int64_t n, GEQRF_MT_TYPE *a, aocl_int64_t lda,
                               GEQRF_MT_TYPE *tau, GEQRF_MT_TYPE *T_storage, uint8_t *a_dependency,
                               aocl_int64_t nb, aocl_int64_t k, GEQRF_MT_TYPE *per_thread_work,
                               aocl_int64_t per_thread_work_size, aocl_int64_t *info)
{
    /* Get the current thread number */
    aocl_int64_t thread_num = omp_get_thread_num();
    aocl_int64_t iinfo;
    /* Get the thread-local workspace */
    GEQRF_MT_TYPE *t_work = per_thread_work + thread_num * per_thread_work_size;

    aocl_int64_t ldwork = nb;

    /* Apply block reflector to all panels in left
     * columns of the current block.
     */
    for(aocl_int64_t i = 0; i < fla_min(m, k); i += nb)
    {
        /* If the block reflector is not ready, wait for it */
        uint8_t ready = 0;
        while(!ready)
        {
#pragma omp atomic read
            ready = a_dependency[i / nb];
            if(!ready)
                /* Use pause instruction for busy waiting */
                __asm__ __volatile__("pause");
        }

        /* Apply the block reflector to the current panel */
        aocl_int64_t mvai = m - i;
        aocl_int64_t n_reflector = fla_min(mvai, nb);
        GEQRF_MT_LARFB_FN("Left", GEQRF_MT_TRANS, "Forward", "Columnwise", &mvai, &n, &n_reflector,
                          &a[i + i * lda], &lda, T_storage + (nb * i), &nb, &a[k * lda + i], &lda,
                          t_work, &ldwork);
    }

    aocl_int64_t mvak = m - k;

    /* If there are no more rows to process, return */
    if(mvak <= 0)
    {
        return;
    }

    /* Factor the current panel */
    GEQRF_MT_GEQR2_FN(&mvak, &n, &a[k * lda + k], &lda, &tau[k], t_work, &iinfo);

    /* If the factorization failed, set the info and return */
    if(iinfo != 0)
    {
        *info = iinfo;
        return;
    }

    /* Generate the triangular factor of the block reflector */
    aocl_int64_t k_reflector = fla_min(mvak, n);
    GEQRF_MT_LARFT_FN("Forward", "Columnwise", &mvak, &k_reflector, &a[k * lda + k], &lda, &tau[k],
                      T_storage + (nb * k), &nb);

    /* Set the block reflector as ready */
#pragma omp atomic write
    a_dependency[k / nb] = 1;
}

static void GEQRF_MT_LARGE_FN(aocl_int64_t gm, aocl_int64_t gn, GEQRF_MT_TYPE *a, aocl_int64_t lda,
                              GEQRF_MT_TYPE *tau, GEQRF_MT_TYPE *work, aocl_int64_t nthreads,
                              aocl_int64_t *info)
{

    aocl_int64_t nb = GEQRF_MT_PANEL_SIZE;

    aocl_int64_t iinfo = 0;

    /* Number of panels */
    aocl_int64_t nt = ((gn - 1) / nb) + 1;

    /* Storage to keep triangular factors for each panel */
    GEQRF_MT_TYPE *T_storage = work;

    /* Dependency list to keep track of which panels are ready */
    uint8_t *a_dependency = (uint8_t *)(T_storage + (nb * nb) * nt);

    /* Per-thread workspace */
    GEQRF_MT_TYPE *per_thread_work
        = ((GEQRF_MT_TYPE *)a_dependency)
          + containers_needed_for_size(nt, sizeof(uint8_t), sizeof(GEQRF_MT_TYPE));

    /* Initialize the dependency list to 0 */
    memset(a_dependency, 0, nt);

    /* Work size for each thread */
    aocl_int64_t per_thread_work_size = (nb * nb);

#if FLA_ENABLE_AOCL_BLAS
    aocl_int64_t orig_blis_threads = bli_thread_get_num_threads();
    bli_thread_set_num_threads(1);
#endif

/* Schedule threads to process panels
 * Panels are assigned to threads in a round-robin manner
 */
#pragma omp parallel for num_threads(nthreads) schedule(static, 1) proc_bind(close)
    for(aocl_int64_t i = 0; i < gn; i += nb)
    {
        /* Get the number of columns in current panel */
        aocl_int64_t n_thread = fla_min(nb, gn - i);
        aocl_int64_t local_info = 0;
        GEQRF_MT_THREAD_FN(gm, n_thread, a, lda, tau, T_storage, a_dependency, nb, i,
                           per_thread_work, per_thread_work_size, &local_info);
        if(local_info != 0)
        {
#pragma omp atomic write
            iinfo = local_info;
        }
    }

#if FLA_ENABLE_AOCL_BLAS
    /* reset no. of threads back to original for BLIS */
    bli_thread_set_num_threads(orig_blis_threads);
#endif

    /* Set the info */
    *info = iinfo;
}

#undef GEQRF_MT_TYPE
#undef GEQRF_MT_GEQR2_FN
#undef GEQRF_MT_LARFB_FN
#undef GEQRF_MT_LARFT_FN
#undef GEQRF_MT_LARGE_FN
#undef GEQRF_MT_THREAD_FN
#undef GEQRF_MT_NUM_THREADS_FN
#undef GEQRF_MT_LWORK_FN
#undef GEQRF_MT_PANEL_SIZE
#undef GEQRF_MT_THRESHOLD_8_THREADS
#undef GEQRF_MT_TRANS
#undef GEQRF_MT_PRECISION

#endif /* FLA_ENABLE_AMD_OPT && FLA_OPENMP_MULTITHREADING */
