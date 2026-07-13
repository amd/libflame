/*
    Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#ifdef FLA_OPENMP_MULTITHREADING

#ifndef API_ID_DEFINED
#define API_ID_DEFINED

#include "FLA_type_defs.h"
#include "FLA_macro_defs.h"
#include <immintrin.h>

/* API ID */
typedef enum
{
    FLA_LABRD = 0,
    FLA_ORMQR,
    FLA_ORMLQ
} API_ID;
#endif

void FLA_Thread_get_subrange(int thread_ID, int num_threads, fla_dim_t range, fla_dim_t *sub_range,
                             fla_dim_t *index);
void FLA_Thread_get_subrange_chunks(int thread_ID, int num_threads,
                                      size_t elem_size_in_bytes, fla_dim_t range, fla_dim_t *sub_range,
                                      fla_dim_t *index, fla_dim_t *thread_threshold);
void FLA_Thread_optimum(API_ID family, int *actual_num_threads);

/*
 * Custom barrier implementation.
 *
 * Notes:
 * - This is a reusable spin barrier intended for tightly-coupled thread teams.
 * - Progress relies on all participants eventually reaching the barrier.
 */

#ifndef FLA_BARRIER_DEFINED
#define FLA_BARRIER_DEFINED

typedef struct __attribute__((aligned(FLA_CACHE_LINE_SIZE_BYTES))) __FLA_BARRIER
{
    /*
     * The two state groups are intentionally placed on separate cache lines.
     * Threads heavily update/poll wait_count during arrival/departure, while
     * release_flag is polled in spin loops. Splitting them reduces false
     * sharing and cache-line ping-pong between counter updates and flag polling.
     */
    /*
     * Arrival state:
     * - total_threads: number of participating threads.
     * - wait_count: number of threads currently inside the barrier.
     */
    struct __attribute__((aligned(FLA_CACHE_LINE_SIZE_BYTES)))
    {
        volatile aocl_int64_t total_threads;
        volatile aocl_int64_t wait_count;
        /* Padding to ensure this struct occupies full cache line */
        char _pad1[FLA_CACHE_LINE_SIZE_BYTES - 2 * sizeof(aocl_int64_t)];
    };
    /*
     * Release state:
     * - release_flag: set to 1 once all threads arrive; reset to 0 by last
     *   thread leaving so the barrier can be reused.
     * - sense: alternating epoch counter to avoid race conditions
     */
    struct __attribute__((aligned(FLA_CACHE_LINE_SIZE_BYTES)))
    {
        volatile aocl_int64_t release_flag;
        volatile aocl_int64_t sense;
        /* Padding to ensure this struct occupies full cache line */
        char _pad2[FLA_CACHE_LINE_SIZE_BYTES - 2 * sizeof(aocl_int64_t)];
    };
} FLA_BARRIER;

/*
 * Initializes a reusable barrier for num_threads participants.
 *
 * Precondition:
 * - num_threads > 0.
 *
 * Typical usage:
 * - Call once before the worker team starts using the barrier.
 * - Do not reinitialize while other threads may still access the barrier.
 */
#define FLA_BARRIER_INIT(barrier, num_threads)                                        \
    do                                                                                \
    {                                                                                 \
        (barrier).total_threads = num_threads; /* Number of participating threads. */ \
        (barrier).wait_count = 0; /* No thread has arrived yet. */                    \
        (barrier).release_flag = 0; /* Barrier starts in closed state. */             \
        (barrier).sense = 0; /* Initial sense/epoch. */                               \
    } while(0)

/*
 * Waits until all participating threads reach the barrier.
 *
 * Sense-reversing barrier (reusable across epochs):
 * 1) Sync: each thread snapshots sense and spins until release_flag equals
 *    that epoch, ensuring the previous round has completed.
 * 2) Arrival: each thread atomically increments wait_count; the last arriver
 *    resets wait_count, bumps sense, and publishes the new epoch in both
 *    sense and release_flag.
 * 3) Release: all other threads spin until release_flag matches the next
 *    epoch (local_sense + 1), then proceed. No separate departure phase.
 */
#define FLA_BARRIER_WAIT(barrier)                                                               \
    do                                                                                          \
    {                                                                                           \
        volatile aocl_int64_t tmp_wait_count; /* Local snapshot of wait_count ops. */           \
        volatile aocl_int64_t local_sense = __atomic_load_n(&(barrier).sense, __ATOMIC_ACQUIRE);\
        int spin_count = 0; /* For exponential backoff. */                                      \
        /* Wait for previous epoch to complete (if any thread still departing). */              \
        while(__atomic_load_n(&(barrier).release_flag, __ATOMIC_ACQUIRE) != local_sense)        \
        {                                                                                       \
            _mm_pause(); /* Reduce contention on release_flag cache line. */                    \
            if(++spin_count > 64)                                                               \
            {                                                                                   \
                spin_count = 0;                                                                 \
                __builtin_ia32_pause();                                                         \
            }                                                                                   \
        }                                                                                       \
        /* Arrival phase: increment counter. */                                                 \
        tmp_wait_count = __atomic_fetch_add(&(barrier).wait_count, 1, __ATOMIC_ACQ_REL);        \
        /* If this is the last arrival, open the gate for all waiters. */                       \
        if(tmp_wait_count == (barrier).total_threads - 1)                                       \
        {                                                                                       \
            /* Reset counter for next epoch. */                                                 \
            __atomic_store_n(&(barrier).wait_count, 0, __ATOMIC_RELEASE);                       \
            /* Flip sense and release all waiting threads. */                                   \
            aocl_int64_t next_sense = local_sense + 1;                                          \
            __atomic_store_n(&(barrier).sense, next_sense, __ATOMIC_RELEASE);                   \
            __atomic_store_n(&(barrier).release_flag, next_sense, __ATOMIC_RELEASE);            \
        }                                                                                       \
        else                                                                                    \
        {                                                                                       \
            /* Wait for last thread to release barrier. */                                      \
            spin_count = 0;                                                                     \
            aocl_int64_t expected_sense = local_sense + 1;                                      \
            while(__atomic_load_n(&(barrier).release_flag, __ATOMIC_ACQUIRE) != expected_sense) \
            {                                                                                   \
                _mm_pause(); /* Spin with pause to reduce memory traffic. */                    \
                if(++spin_count > 64)                                                           \
                {                                                                               \
                    spin_count = 0;                                                             \
                    for(int i = 0; i < 8; i++)                                                  \
                        _mm_pause();                                                            \
                }                                                                               \
            }                                                                                   \
        }                                                                                       \
    } while(0)

#endif /* FLA_BARRIER_DEFINED */

/*
 * Atomically performs: *ptr = max(*ptr, max_val).
 *
 * Algorithm:
 * - Load current value.
 * - If current value is already >= max_val, return.
 * - Otherwise attempt CAS to publish max_val.
 * - On CAS failure, retry using the refreshed value returned via old_val.
 *
 * Memory-order intent:
 * - RELEASE on successful CAS publishes prior writes before the new maximum.
 * - RELAXED on load/failure path is sufficient for retry mechanics.
 */
#define fla_atomic_max(ptr, max_val)                                                             \
    do                                                                                           \
    {                                                                                            \
        aocl_int64_t old_val = __atomic_load_n((ptr), __ATOMIC_RELAXED); /* Initial snapshot. */ \
        /* CAS succeeds only if *ptr is still old_val; otherwise old_val is refreshed. */        \
        while(old_val < (max_val)                                                                \
              && !__atomic_compare_exchange_n((ptr), &old_val, (max_val), 0, __ATOMIC_RELEASE,   \
                                              __ATOMIC_RELAXED))                                 \
            ;                                                                                    \
    } while(0)

#endif