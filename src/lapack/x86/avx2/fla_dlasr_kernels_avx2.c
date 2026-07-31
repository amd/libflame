/******************************************************************************
 * * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* The side='L' rotation loops below walk columns with stride a_dim1 while
   touching adjacent rows, which are contiguous. Left to itself GCC 15+ merges
   the adjacent-row accesses into contiguous vector accesses and pays for an
   in-register transpose of cross-lane shuffles, which is throughput-bound on
   Zen and costs ~3x. Assembling each column pair explicitly with a low/high
   load pair keeps the data in lanes and leaves the vectorizer nothing to
   re-vectorize.

   A 256-bit ("4 columns per iteration") variant that assembles its vectors
   the same way (no gather) was benchmarked and does give a real 15-45%
   speedup whenever the per-rotation working set (n columns, one cache line
   each) is small enough to stay resident in L1/L2 -- i.e. it helps for the
   common small/medium n case and is a wash once n is large enough to be
   memory-bandwidth-bound. AVX512 was also tried (gather/scatter and a
   512-bit insert variant) and gave no benefit at all -- worse with gather,
   a wash with insert -- so there is no 512-bit kernel here.

   One caveat: whenever a_dim1 (lda) is an exact multiple of 64 elements,
   *both* kernels hit a severe, pre-existing cache/set-associativity cliff
   (10-40x slower than neighboring lda values -- classic "power-of-2 leading
   dimension" pathology, independent of this file) and the 256- vs 128-bit
   comparison becomes erratic in that narrow band (sometimes the wider one
   wins bigger, sometimes it loses by ~20%). fla_dlasr_left_pivotv_avx2()
   below simply keeps to the 128-bit kernel whenever a_dim1 is a multiple of
   64, since that is already what ships today and avoids adding new
   variance to an already pathological case. */

/* Apply one plane rotation (ct,st) to the adjacent row pair at base[0],base[1]
   across n columns. */
static inline void fla_dlasr_rot1_avx2(doublereal *base, aocl_int64_t n, aocl_int64_t a_dim1,
                                       doublereal ct, doublereal st)
{
    const __m128d vct = _mm_set1_pd(ct);
    const __m128d vst = _mm_set1_pd(st);
    aocl_int64_t i = 1;

    for(; i + 1 <= n; i += 2)
    {
        doublereal *q0 = base + i * a_dim1;
        doublereal *q1 = q0 + a_dim1;
        __m128d t0 = _mm_loadh_pd(_mm_load_sd(q0), q1);
        __m128d t1 = _mm_loadh_pd(_mm_load_sd(q0 + 1), q1 + 1);
        __m128d r0, r1;

        r1 = _mm_fmsub_pd(vct, t1, _mm_mul_pd(vst, t0));
        r0 = _mm_fmadd_pd(vst, t1, _mm_mul_pd(vct, t0));

        _mm_storel_pd(q0, r0);
        _mm_storeh_pd(q1, r0);
        _mm_storel_pd(q0 + 1, r1);
        _mm_storeh_pd(q1 + 1, r1);
    }

    for(; i <= n; ++i)
    {
        doublereal *q0 = base + i * a_dim1;
        doublereal temp = q0[1];
        q0[1] = ct * temp - st * q0[0];
        q0[0] = st * temp + ct * q0[0];
    }
}

/* Apply two consecutive plane rotations to the three adjacent rows at base[0],
   base[1], base[2] across n columns. The second rotation consumes the row
   updated by the first, so the pair is fused to keep that value in register. */
static inline void fla_dlasr_rot2_avx2(doublereal *base, aocl_int64_t n, aocl_int64_t a_dim1,
                                       doublereal ct0, doublereal st0, doublereal ct1,
                                       doublereal st1)
{
    const __m128d vct0 = _mm_set1_pd(ct0);
    const __m128d vst0 = _mm_set1_pd(st0);
    const __m128d vct1 = _mm_set1_pd(ct1);
    const __m128d vst1 = _mm_set1_pd(st1);
    aocl_int64_t i = 1;

    for(; i + 1 <= n; i += 2)
    {
        doublereal *q0 = base + i * a_dim1;
        doublereal *q1 = q0 + a_dim1;
        __m128d t0 = _mm_loadh_pd(_mm_load_sd(q0), q1);
        __m128d t1 = _mm_loadh_pd(_mm_load_sd(q0 + 1), q1 + 1);
        __m128d t2 = _mm_loadh_pd(_mm_load_sd(q0 + 2), q1 + 2);
        __m128d r0, r1, r2, m1;

        r0 = _mm_fmadd_pd(vst0, t1, _mm_mul_pd(vct0, t0));
        m1 = _mm_fmsub_pd(vct0, t1, _mm_mul_pd(vst0, t0));
        r1 = _mm_fmadd_pd(vst1, t2, _mm_mul_pd(vct1, m1));
        r2 = _mm_fmsub_pd(vct1, t2, _mm_mul_pd(vst1, m1));

        _mm_storel_pd(q0, r0);
        _mm_storeh_pd(q1, r0);
        _mm_storel_pd(q0 + 1, r1);
        _mm_storeh_pd(q1 + 1, r1);
        _mm_storel_pd(q0 + 2, r2);
        _mm_storeh_pd(q1 + 2, r2);
    }

    for(; i <= n; ++i)
    {
        doublereal *q0 = base + i * a_dim1;
        doublereal tmp0 = q0[0], tmp1 = q0[1], tmp2 = q0[2];
        doublereal res0 = ct0 * tmp0 + st0 * tmp1;

        tmp1 = ct0 * tmp1 - st0 * tmp0;
        q0[0] = res0;
        q0[1] = ct1 * tmp1 + st1 * tmp2;
        q0[2] = ct1 * tmp2 - st1 * tmp1;
    }
}

/* 256-bit counterpart of fla_dlasr_rot1_avx2: same idea, 4 columns per
 * iteration. Any trailing columns (n not a multiple of 4) are handed off to
 * the 128-bit kernel above, which already knows how to close out with its
 * own 2-wide/scalar tail. */
static inline void fla_dlasr_rot1_avx2_256(doublereal *base, aocl_int64_t n, aocl_int64_t a_dim1,
                                           doublereal ct, doublereal st)
{
    const __m256d vct = _mm256_set1_pd(ct);
    const __m256d vst = _mm256_set1_pd(st);
    aocl_int64_t i = 1;

    for(; i + 3 <= n; i += 4)
    {
        doublereal *q0 = base + i * a_dim1;
        doublereal *q1 = q0 + a_dim1;
        doublereal *q2 = q1 + a_dim1;
        doublereal *q3 = q2 + a_dim1;

        __m128d p0 = _mm_loadu_pd(q0);
        __m128d p1 = _mm_loadu_pd(q1);
        __m128d p2 = _mm_loadu_pd(q2);
        __m128d p3 = _mm_loadu_pd(q3);

        __m128d lo01 = _mm_unpacklo_pd(p0, p1);
        __m128d hi01 = _mm_unpackhi_pd(p0, p1);
        __m128d lo23 = _mm_unpacklo_pd(p2, p3);
        __m128d hi23 = _mm_unpackhi_pd(p2, p3);

        __m256d t0 = _mm256_set_m128d(lo23, lo01);
        __m256d t1 = _mm256_set_m128d(hi23, hi01);
        __m256d r0, r1;

        r1 = _mm256_fmsub_pd(vct, t1, _mm256_mul_pd(vst, t0));
        r0 = _mm256_fmadd_pd(vst, t1, _mm256_mul_pd(vct, t0));

        __m128d r0_lo = _mm256_castpd256_pd128(r0);
        __m128d r0_hi = _mm256_extractf128_pd(r0, 1);
        __m128d r1_lo = _mm256_castpd256_pd128(r1);
        __m128d r1_hi = _mm256_extractf128_pd(r1, 1);

        _mm_storeu_pd(q0, _mm_unpacklo_pd(r0_lo, r1_lo));
        _mm_storeu_pd(q1, _mm_unpackhi_pd(r0_lo, r1_lo));
        _mm_storeu_pd(q2, _mm_unpacklo_pd(r0_hi, r1_hi));
        _mm_storeu_pd(q3, _mm_unpackhi_pd(r0_hi, r1_hi));
    }

    if(i <= n)
    {
        fla_dlasr_rot1_avx2(base + (i - 1) * a_dim1, n - i + 1, a_dim1, ct, st);
    }
}

/* 256-bit counterpart of fla_dlasr_rot2_avx2: same fused pair of rotations,
 * 4 columns per iteration. The third row (base[2]) only contributes one
 * double per column, so it is assembled/stored as a lane vector rather than
 * via the unpack trick used for the row pairs. */
static inline void fla_dlasr_rot2_avx2_256(doublereal *base, aocl_int64_t n, aocl_int64_t a_dim1,
                                           doublereal ct0, doublereal st0, doublereal ct1,
                                           doublereal st1)
{
    const __m256d vct0 = _mm256_set1_pd(ct0);
    const __m256d vst0 = _mm256_set1_pd(st0);
    const __m256d vct1 = _mm256_set1_pd(ct1);
    const __m256d vst1 = _mm256_set1_pd(st1);
    aocl_int64_t i = 1;

    for(; i + 3 <= n; i += 4)
    {
        doublereal *q0 = base + i * a_dim1;
        doublereal *q1 = q0 + a_dim1;
        doublereal *q2 = q1 + a_dim1;
        doublereal *q3 = q2 + a_dim1;

        __m128d p0 = _mm_loadu_pd(q0);
        __m128d p1 = _mm_loadu_pd(q1);
        __m128d p2 = _mm_loadu_pd(q2);
        __m128d p3 = _mm_loadu_pd(q3);

        __m128d lo01 = _mm_unpacklo_pd(p0, p1);
        __m128d hi01 = _mm_unpackhi_pd(p0, p1);
        __m128d lo23 = _mm_unpacklo_pd(p2, p3);
        __m128d hi23 = _mm_unpackhi_pd(p2, p3);

        __m256d t0 = _mm256_set_m128d(lo23, lo01);
        __m256d t1 = _mm256_set_m128d(hi23, hi01);
        __m256d t2 = _mm256_set_pd(q3[2], q2[2], q1[2], q0[2]);
        __m256d r0, r1, r2, m1;

        r0 = _mm256_fmadd_pd(vst0, t1, _mm256_mul_pd(vct0, t0));
        m1 = _mm256_fmsub_pd(vct0, t1, _mm256_mul_pd(vst0, t0));
        r1 = _mm256_fmadd_pd(vst1, t2, _mm256_mul_pd(vct1, m1));
        r2 = _mm256_fmsub_pd(vct1, t2, _mm256_mul_pd(vst1, m1));

        __m128d r0_lo = _mm256_castpd256_pd128(r0);
        __m128d r0_hi = _mm256_extractf128_pd(r0, 1);
        __m128d r1_lo = _mm256_castpd256_pd128(r1);
        __m128d r1_hi = _mm256_extractf128_pd(r1, 1);

        _mm_storeu_pd(q0, _mm_unpacklo_pd(r0_lo, r1_lo));
        _mm_storeu_pd(q1, _mm_unpackhi_pd(r0_lo, r1_lo));
        _mm_storeu_pd(q2, _mm_unpacklo_pd(r0_hi, r1_hi));
        _mm_storeu_pd(q3, _mm_unpackhi_pd(r0_hi, r1_hi));

        doublereal r2v[4];
        _mm256_storeu_pd(r2v, r2);
        q0[2] = r2v[0];
        q1[2] = r2v[1];
        q2[2] = r2v[2];
        q3[2] = r2v[3];
    }

    if(i <= n)
    {
        fla_dlasr_rot2_avx2(base + (i - 1) * a_dim1, n - i + 1, a_dim1, ct0, st0, ct1, st1);
    }
}

/* Apply a sequence of plane rotations from the left with a variable pivot
 * (side='L', pivot='V'); forward != 0 selects direct='F', else direct='B'.
 * c__, s and a are the f2c-adjusted one-based pointers, so the callers'
 * --c__, --s and a -= a_offset must already have been applied.
 * */
void fla_dlasr_left_pivotv_avx2(logical forward, aocl_int64_t m, aocl_int64_t n, doublereal *c__,
                                doublereal *s, doublereal *a, aocl_int64_t a_dim1)
{
    aocl_int64_t j, mm1;
    doublereal ctemp, stemp;
    /* See the comment at the top of this file: lda that is an exact
       multiple of 64 elements is a pre-existing cache-associativity cliff
       for both kernel widths, and the two widths behave erratically
       relative to each other in that band. Keep to the 128-bit kernel
       there since it is what already ships. */
    logical use_256 = ((a_dim1 & 63) != 0);

    if(forward)
    {
        mm1 = m - 1;
        /* Apply two rotations in an iteration */
        for(j = 1; j < mm1; j += 2)
        {
            if(use_256)
            {
                fla_dlasr_rot2_avx2_256(a + j, n, a_dim1, c__[j], s[j], c__[j + 1], s[j + 1]);
            }
            else
            {
                fla_dlasr_rot2_avx2(a + j, n, a_dim1, c__[j], s[j], c__[j + 1], s[j + 1]);
            }
        }
        /* Apply the remaining rotation */
        if(mm1 & 1)
        {
            ctemp = c__[mm1];
            stemp = s[mm1];
            if(ctemp != 1. || stemp != 0.)
            {
                if(use_256)
                {
                    fla_dlasr_rot1_avx2_256(a + mm1, n, a_dim1, ctemp, stemp);
                }
                else
                {
                    fla_dlasr_rot1_avx2(a + mm1, n, a_dim1, ctemp, stemp);
                }
            }
        }
    }
    else
    {
        for(j = m - 1; j >= 1; --j)
        {
            ctemp = c__[j];
            stemp = s[j];
            if(ctemp != 1. || stemp != 0.)
            {
                if(use_256)
                {
                    fla_dlasr_rot1_avx2_256(a + j, n, a_dim1, ctemp, stemp);
                }
                else
                {
                    fla_dlasr_rot1_avx2(a + j, n, a_dim1, ctemp, stemp);
                }
            }
        }
    }
    return;
}
#endif
