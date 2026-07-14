/*
   Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*
 * Tile-based, lock-free parallel Cholesky factorization (DPOTRF).
 *
 * The matrix is split into NB x NB tiles. Each tile has an atomic "step"
 * counter tracking how far it has been processed. Worker threads scan the
 * grid and claim any tile whose dependencies are met. No locks, no queue:
 * ordering comes entirely from the per-tile counters and acquire/release
 * atomics.
 *
 * Tile dependencies (uplo='L', A = L * L^T). For diagonal tile i, every tile
 * in panel i must absorb a SYRK/GEMM update from each earlier panel before it
 * can be finalized:
 *
 *   [L00]
 *   [L10] [L11]                            L11 = potrf(A11 - L10*L10^T)
 *   [L20] [L21] [L22]                      L22 = potrf(A22 - L20*L20^T - L21*L21^T)
 *   [L30] [L31] [L32] [L33]                L33 = potrf(A33 - ... )
 *
 *   panel 0:  potrf(A00) -> L00,   then trsm to get L10,L20,L30
 *   panel 1:  syrk/gemm updates from panel 0,  potrf -> L11,  trsm -> L21,L31
 *   panel 2:  updates from panels 0,1,         potrf -> L22,  trsm -> L32
 *
 * Per-tile `step` counter lifecycle. For a tile in panel i, `step` counts how
 * many updates it has absorbed; it is finalized at step == i and then marked
 * done as i+1:
 *
 *   step:  0      1      ...    i              i+1
 *          |      |             |               |
 *          +--update from-------+--finalize-----+--> done
 *           panels 0..i-1        (potrf if i==j,
 *           (syrk/gemm)           else trsm)
 *
 * The BUSY bit is OR'd into `step` while a thread owns the tile, so exactly
 * one thread advances it per step.
 */

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#ifdef FLA_OPENMP_MULTITHREADING
#include <omp.h>
#endif
#if FLA_ENABLE_AMD_OPT
int fla_dpotrf_small_avx2(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                          aocl_int64_t *info);
#endif

/* Minimum matrix order to use the parallel path. */
#define FLA_DPOTRF_THRESHOLD 256
/* Spacing between counters so each sits on its own cache line (no false sharing). */
#define FLA_DPOTRF_PAD 8
/* Bit set in a step counter to mark a tile as claimed. */
#define FLA_DPOTRF_BUSY ((aocl_int64_t)1 << 62)
extern int lapack_dpotrf(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                         aocl_int64_t *info);
extern int fla_thread_get_num_threads(void);

/* BLAS scalars passed by address. */
doublereal one = 1.0, neg_one = -1.0;
/* Shared context, set up once then read by every worker thread. */
typedef struct
{
    doublereal *a; /* matrix base pointer                     */
    aocl_int64_t lda; /* leading dimension                       */
    aocl_int64_t n; /* matrix order                            */
    aocl_int64_t nb; /* tile block size                         */
    aocl_int64_t nt; /* tiles per dimension                     */
    aocl_int64_t total; /* tiles in the active triangle            */
    int upper; /* 1 for 'U', 0 for 'L'                    */
    aocl_int64_t *step; /* per-tile step counters (padded)         */
    aocl_int64_t *coldone; /* finalized-tile count per panel (padded) */
    aocl_int64_t done; /* atomic count of completed tiles         */
    aocl_int64_t err; /* atomic first failure position (>0)      */
} dpotrf_ctx;

/* Pick (num_threads, block_size) per size band. */
static void dpotrf_auto_tune_params(aocl_int64_t n, int *num_threads, aocl_int64_t *block_size)
{
    int max_threads = fla_thread_get_num_threads();
    int want;
    aocl_int64_t nb;

    if(n <= 300)
    {
        want = 8;
        nb = 64;
    }
    else if(n <= 1024)
    {
        want = 8;
        nb = 64;
    }
    else if(n < 3000)
    {
        want = 32;
        nb = 96;
    }
    else if(n <= 4800)
    {
        want = 64;
        nb = 128;
    }
    else if(n < 8192)
    {
        want = 128;
        nb = 192;
    }
    else if(n < 9000)
    {
        want = 192;
        nb = 192;
    }
    else
    {
        want = 128;
        nb = 192;
    }

    if(want > max_threads)
        want = max_threads;
    if(want < 1)
        want = 1;

    *num_threads = want;
    *block_size = nb;
}
/* Size of a tile. NB or the remainder for the last tile. */
#define TILE_DIM(ctx, idx) \
    (((idx) + 1) * (ctx)->nb <= (ctx)->n ? (ctx)->nb : (ctx)->n - (idx) * (ctx)->nb)
/* Address of a tile (row, col). */
#define TILE_PTR(ctx, row, col) (&(ctx)->a[(row) * (ctx)->nb + (col) * (ctx)->nb * (ctx)->lda])

/* Map a (i, j) tile to its index in the padded `step` array.
 * Both uplo cases are folded onto the same triangular layout. */
static inline aocl_int64_t get_step_index(dpotrf_ctx *ctx, aocl_int64_t i, aocl_int64_t j)
{
    return (ctx->upper ? (i * ctx->nt + j) : (j * ctx->nt + i)) * FLA_DPOTRF_PAD;
}

/* Index into the padded `coldone` array for panel i. */
static inline aocl_int64_t get_coldone_index(aocl_int64_t i)
{
    return i * FLA_DPOTRF_PAD;
}

/* Factorize a diagonal tile (tuned AVX2 kernel when available). */
inline static void potrf_factor(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                                aocl_int64_t *info)
{
#if FLA_ENABLE_AMD_OPT
    fla_dpotrf_small_avx2(uplo, n, a, lda, info);
#else
    aocl_lapack_dpotrf(uplo, n, a, lda, info);
#endif
}

/*
 * Run one tile operation on (i, j). `step` selects which one:
 *   step <  i : trailing update from panel `step` (SYRK if i==j, else GEMM).
 *   step == i : finalize the tile (POTRF if i==j, else TRSM).
 * A tile takes i updates before being finalized at step==i.
 */
static void dpotrf_execute(dpotrf_ctx *ctx, aocl_int64_t i, aocl_int64_t j, aocl_int64_t step)
{
    aocl_int64_t lda = ctx->lda;
    aocl_int64_t ldb = ctx->lda;
    aocl_int64_t ldc = ctx->lda;
    if(step < i)
    {
        /* Trailing update: subtract the contribution of panel `step`. */
        aocl_int64_t k = TILE_DIM(ctx, step);
        if(i == j)
        {
            // syrk: update the diagonal tile
            aocl_int64_t n = TILE_DIM(ctx, i);
            doublereal *C = TILE_PTR(ctx, i, i);
            if(ctx->upper)
            {
                doublereal *A = TILE_PTR(ctx, step, i);
                aocl_blas_dsyrk("Upper", "Transpose", &n, &k, &neg_one, A, &lda, &one, C, &ldc);
            }
            else
            {
                doublereal *A = TILE_PTR(ctx, i, step);
                aocl_blas_dsyrk("Lower", "No Transpose", &n, &k, &neg_one, A, &lda, &one, C, &ldc);
            }
        }
        else
        {
            // gemm: update an off-diagonal tile
            aocl_int64_t m = TILE_DIM(ctx, i);
            aocl_int64_t n = TILE_DIM(ctx, j);
            if(ctx->upper)
            {
                doublereal *A = TILE_PTR(ctx, step, i);
                doublereal *B = TILE_PTR(ctx, step, j);
                doublereal *C = TILE_PTR(ctx, i, j);
                aocl_blas_dgemm("Transpose", "No Transpose", &m, &n, &k, &neg_one, A, &lda, B, &ldb,
                                &one, C, &ldc);
            }
            else
            {
                doublereal *A = TILE_PTR(ctx, j, step);
                doublereal *B = TILE_PTR(ctx, i, step);
                doublereal *C = TILE_PTR(ctx, j, i);
                aocl_blas_dgemm("No Transpose", "Transpose", &n, &m, &k, &neg_one, A, &lda, B, &ldb,
                                &one, C, &ldc);
            }
        }
        return;
    }
    if(i == j)
    {
        // potrf: factorize the fully-updated diagonal tile
        char uplo = ctx->upper ? 'U' : 'L';
        aocl_int64_t n = TILE_DIM(ctx, i);
        doublereal *A = TILE_PTR(ctx, i, i);
        aocl_int64_t info_ = 0;
        potrf_factor(&uplo, &n, A, &lda, &info_);
        if(info_ != 0)
        {
            /* Not positive definite: record the first error position */
            aocl_int64_t expected = 0;
            aocl_int64_t position_of_error = i * ctx->nb + info_;
            __atomic_compare_exchange_n(&ctx->err, &expected, position_of_error, 0,
                                        __ATOMIC_ACQ_REL, __ATOMIC_RELAXED);
        }
    }
    else
    {
        /* trsm: apply the factored diagonal block to its panel tile */
        doublereal *A = TILE_PTR(ctx, i, i);
        aocl_int64_t m = TILE_DIM(ctx, i);
        aocl_int64_t n = TILE_DIM(ctx, j);
        if(ctx->upper)
        {
            /* A(i,j) = U(i,i)^{-T} * A(i,j) */
            doublereal *B = TILE_PTR(ctx, i, j);
            aocl_blas_dtrsm("Left", "Upper", "Transpose", "Non-unit", &m, &n, &one, A, &lda, B,
                            &ldb);
        }
        else
        {
            /* A(j,i) = A(j,i) * L(i,i)^{-T} */
            doublereal *B = TILE_PTR(ctx, j, i);
            aocl_blas_dtrsm("Right", "Lower", "Transpose", "Non-unit", &n, &m, &one, A, &lda, B,
                            &ldb);
        }
    }
}

/* True once the producer tile (step, i) has advanced to step+1, i.e. the
 * data this consumer needs is ready and visible. */
inline static int dpotrf_producer_done(dpotrf_ctx *ctx, aocl_int64_t i, aocl_int64_t step)
{
    aocl_int64_t completed_step
        = __atomic_load_n(&ctx->step[get_step_index(ctx, step, i)], __ATOMIC_ACQUIRE);
    return completed_step == step + 1;
}

/*
 * Try to advance tile (i, j) by one step. Returns 1 if this thread did the
 * work, 0 if the tile wasn't ready, was done, or another thread claimed it.
 * Steps: check ready -> CAS the BUSY bit to claim -> run -> publish next step.
 */
static int dpotrf_try(dpotrf_ctx *ctx, aocl_int64_t i, aocl_int64_t j)
{
    aocl_int64_t index = get_step_index(ctx, i, j);
    aocl_int64_t step = __atomic_load_n(&ctx->step[index], __ATOMIC_ACQUIRE);
    aocl_int64_t ready = 0;
    /* Skip if claimed by another thread or already finalized. */
    if((step & FLA_DPOTRF_BUSY) || step == i + 1)
        return 0;
    if(step < i)
    {
        /* Trailing update: needs one producer (diagonal) or two (off-diagonal). */
        if(i == j)
            ready = dpotrf_producer_done(ctx, i, step);
        else
            ready = dpotrf_producer_done(ctx, i, step) && dpotrf_producer_done(ctx, j, step);
    }
    else
    {
        /* Finalization: POTRF is always ready; TRSM needs its diagonal factored. */
        ready = (j == i) ? 1
                         : (__atomic_load_n(&ctx->step[get_step_index(ctx, i, i)], __ATOMIC_ACQUIRE)
                            == i + 1);
    }
    if(!ready)
        return 0;
    aocl_int64_t expected = step;
    aocl_int64_t locked = step | FLA_DPOTRF_BUSY;

    /* Claim the tile; only one thread can flip BUSY from `step`. */
    if(!__atomic_compare_exchange_n(&ctx->step[index], &expected, locked, 0, __ATOMIC_ACQ_REL,
                                    __ATOMIC_ACQUIRE))
        return 0;

    dpotrf_execute(ctx, i, j, step);

    /* Publish the next step: i+1 if finalized, else step+1. This also clears
     * BUSY and makes the result visible to other threads. */
    aocl_int64_t next = step + 1;
    __atomic_store_n(&ctx->step[index], next, __ATOMIC_RELEASE);
    if(next == i + 1)
    {
        __atomic_fetch_add(&ctx->done, 1, __ATOMIC_ACQ_REL);
        __atomic_fetch_add(&ctx->coldone[get_coldone_index(i)], 1, __ATOMIC_RELEASE);
    }
    return 1;
}

/*
 * Worker loop run by every thread. Repeatedly scans the grid for ready tiles
 * and claims them via dpotrf_try(): first the critical path (diagonal and its
 * panel), then any trailing update. Each thread starts at a different offset
 * (`rotation`) to reduce contention. Exits on error or full completion.
 */
static void dpotrf_worker(dpotrf_ctx *ctx)
{
    aocl_int64_t nt = ctx->nt;
#ifdef FLA_OPENMP_MULTITHREADING
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    int progressed = 0;
    aocl_int64_t i, ii, j, jj;
    /* Per-thread scan offset to spread work and reduce contention. */
    aocl_int64_t rotation = nt > 0 ? (aocl_int64_t)(tid % nt) : 0;

    /* Lowest panel not yet fully finalized; scans skip everything below it. */
    aocl_int64_t completed_column = 0;

    while(1)
    {
        /* Stop on failure or once all tiles are done. */
        if(__atomic_load_n(&ctx->err, __ATOMIC_ACQUIRE) != 0)
            return;
        if(__atomic_load_n(&ctx->done, __ATOMIC_ACQUIRE) >= ctx->total)
            return;
        progressed = 0;
        /* Skip past panels that are fully finalized. Panel p has (nt - p)
         * tiles, so it's done when coldone[p] reaches that count. */
        while(
            completed_column < nt
            && __atomic_load_n(&ctx->coldone[get_coldone_index(completed_column)], __ATOMIC_ACQUIRE)
                   == nt - completed_column)
            completed_column++;

        /* Phase 1: critical path. Factor each diagonal as soon as it's ready,
         * then immediately do its panel TRSM tiles so consumers unblock early.
         * Stop at the first unfactored diagonal (later panels can't be ahead). */
        for(i = completed_column; i < nt; i++)
        {
            if(dpotrf_try(ctx, i, i))
            {
                progressed = 1;
                break;
            }
            if(__atomic_load_n(&ctx->step[get_step_index(ctx, i, i)], __ATOMIC_ACQUIRE) == i + 1)
            {
                /* Diagonal i is factored: finalize its panel tiles. */
                aocl_int64_t tiles_per_column = nt - (i + 1);
                for(j = 0; j < tiles_per_column; j++)
                {
                    aocl_int64_t j_ = i + 1 + (j + rotation) % tiles_per_column;
                    /* Step == i means this panel tile is ready for its TRSM. */
                    if(__atomic_load_n(&ctx->step[get_step_index(ctx, i, j_)], __ATOMIC_ACQUIRE)
                           == i
                       && dpotrf_try(ctx, i, j_))
                    {
                        progressed = 1;
                        break;
                    }
                }
                if(progressed)
                    break;
            }
            else
            {
                break;
            }
        }
        if(progressed)
            continue;
        /* Phase 2: no critical-path work, so do any pending trailing update
         * (SYRK/GEMM, step < i_) in the triangle, again from the rotation offset. */
        aocl_int64_t tiles_per_row = nt - completed_column;
        for(ii = 0; ii < tiles_per_row && !progressed; ii++)
        {
            aocl_int64_t i_ = completed_column + (ii + rotation) % tiles_per_row;
            for(jj = i_; jj < nt; jj++)
            {
                aocl_int64_t step
                    = __atomic_load_n(&ctx->step[get_step_index(ctx, i_, jj)], __ATOMIC_ACQUIRE);
                if(step >= i_)
                    continue; /* only a finalization left; handled in phase 1 */
                if(dpotrf_try(ctx, i_, jj))
                {
                    progressed = 1;
                    break;
                }
            }
        }
        if(!progressed)
        {
            /* Nothing to do this pass: pause before rescanning. */
            __asm__ __volatile__("pause");
        }
    }
}

/*
 * Set up the context, run the worker team in place on `a`, and return the
 * first error position in *info (0 on success). Falls back to reference DPOTRF
 * if allocations fail.
 */
static int dpotrf_var2(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                       aocl_int64_t *info)
{
    dpotrf_ctx ctx;
    int num_threads;
    ctx.a = a;
    ctx.lda = *lda;
    ctx.n = *n;
    dpotrf_auto_tune_params(ctx.n, &num_threads, &ctx.nb);
    /* Tiles per dimension (rounded up) and tiles in the triangle. */
    ctx.nt = ctx.n / ctx.nb + (ctx.n % ctx.nb != 0);
    ctx.total = ctx.nt * (ctx.nt + 1) / 2;
    ctx.upper = (uplo[0] == 'U' || uplo[0] == 'u') ? 1 : 0;
    /* calloc starts every counter at 0 (= "not started"). */
    ctx.step
        = (aocl_int64_t *)calloc((size_t)ctx.nt * ctx.nt * FLA_DPOTRF_PAD, sizeof(aocl_int64_t));
    ctx.coldone = (aocl_int64_t *)calloc((size_t)ctx.nt * FLA_DPOTRF_PAD, sizeof(aocl_int64_t));
    if(ctx.step == NULL || ctx.coldone == NULL)
    {
        free(ctx.step);
        free(ctx.coldone);
        return lapack_dpotrf(uplo, n, a, lda, info);
    }
    ctx.done = 0;
    ctx.err = 0;

    /* The barrier at the end of the parallel region ensures all work is done
     * before we read ctx.err. Guard the pragma so builds without OpenMP
     * (e.g. multithreading disabled) don't emit -Wunknown-pragmas; the block
     * below is valid with or without the parallel region. */
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel num_threads(num_threads) shared(ctx)
#else
    (void)num_threads;
#endif
    {
        dpotrf_worker(&ctx);
    }

    *info = ctx.err;
    free(ctx.step);
    free(ctx.coldone);
    return 0;
}

/*
 * Public entry point. Validates arguments, then uses the tile-parallel path
 * for large enough multithreaded problems, else the reference blocked DPOTRF.
 */
int lapack_dpotrf_var2(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                       aocl_int64_t *info)
{

    logical upper;
    aocl_int64_t i__1;
#ifndef FLA_ENABLE_AOCL_BLAS
    extern logical lsame_(char *ca, char *cb, aocl_int64_t a, aocl_int64_t b);
#endif
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DPOTRF", &i__1, (ftnlen)6);
        return 0;
    }

    if(*n == 0)
    {
        return 0;
    }

#ifdef FLA_OPENMP_MULTITHREADING
    {
        /* Parallel path only pays off with >1 thread and a large matrix. */
        int num_threads = fla_thread_get_num_threads();
        if(num_threads > 1 && *n >= FLA_DPOTRF_THRESHOLD)
        {
            return dpotrf_var2(uplo, n, a, lda, info);
        }
    }
#endif

    /* Small problems / single-threaded: reference blocked path. */
    return lapack_dpotrf(uplo, n, a, lda, info);
}