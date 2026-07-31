/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_GEQRF_SMALL_AVX2_KERNEL_H
#define FLA_GEQRF_SMALL_AVX2_KERNEL_H

/*! @file fla_geqrf_small_avx2_kernel.h
 *  @brief QR kernels for small sizes (real s/d).
 *  */

#if FLA_ENABLE_AMD_OPT

#include "immintrin.h"

#define FLA_GEQRF_SMALL_CAT_(a, b) a##b
#define FLA_GEQRF_SMALL_CAT(a, b) FLA_GEQRF_SMALL_CAT_(a, b)

#define d_ELT doublereal
#define d_V2 __m128d
#define d_V4 __m256d
#define d_W2 2
#define d_W4 4

#define d_THRES_SML 1.491668e-154
#define d_THRES_BIG 1.997919e+146
#define d_SCALE_SML 4.498914e+161
#define d_SCALE_BIG 1.111379e-162

#define d_LAMCH(s) dlamch_(s)
#define d_SCAL aocl_blas_dscal

#define d_V2_SET1(x) _mm_set1_pd(x)
#define d_V4_SET1(x) _mm256_set1_pd(x)
#define d_V2_SETZERO() _mm_setzero_pd()
#define d_V4_SETZERO() _mm256_setzero_pd()
#define d_V2_LOADU(p) _mm_loadu_pd((const doublereal *)(p))
#define d_V4_LOADU(p) _mm256_loadu_pd((const doublereal *)(p))
#define d_V2_STOREU(p, v) _mm_storeu_pd((doublereal *)(p), v)
#define d_V4_STOREU(p, v) _mm256_storeu_pd((doublereal *)(p), v)
#define d_V2_LOAD_S(p) _mm_load_sd((const doublereal *)(p))
#define d_V2_STORE_S(p, v) _mm_store_sd((doublereal *)(p), v)
#define d_V2_LOAD_DUP(p) _mm_loaddup_pd((const doublereal *)(p))
#define d_V2_STORE_L(p, v) _mm_storel_pd((doublereal *)(p), v)
#define d_V2_LOAD_DUP_NORM(p) _mm_loaddup_pd((const doublereal *)(p))
#define d_V2_MUL(a, b) _mm_mul_pd(a, b)
#define d_V4_MUL(a, b) _mm256_mul_pd(a, b)
#define d_V2_ADD(a, b) _mm_add_pd(a, b)
#define d_V4_ADD(a, b) _mm256_add_pd(a, b)
#define d_V2_ANDNOT(a, b) _mm_andnot_pd(a, b)
#define d_V4_ANDNOT(a, b) _mm256_andnot_pd(a, b)
#define d_V2_CMP(a, b, c) _mm_cmp_pd(a, b, c)
#define d_V4_CMP(a, b, c) _mm256_cmp_pd(a, b, c)
#define d_V2_OR(a, b) _mm_or_pd(a, b)
#define d_V4_OR(a, b) _mm256_or_pd(a, b)
#define d_V2_TESTZ(a, b) _mm_testz_pd(a, b)
#define d_V4_TESTZ(a, b) _mm256_testz_pd(a, b)
#define d_V2_BLENDV(a, b, m) _mm_blendv_pd(a, b, m)
#define d_V4_BLENDV(a, b, m) _mm256_blendv_pd(a, b, m)
#define d_V2_HSUM(v) ((v)[0] + (v)[1])
#define d_V4_HSUM(v) ((v)[0] + (v)[1] + (v)[2] + (v)[3])
#define d_TAIL_LOOP(x, lim) if((x) == (lim))

#define d_INIT_EXTRA __m128d vd2_ltmp, vd2_htmp; \
    __m128d vd2_ntau;
#define d_APPLY_NTAU_SETUP vd2_ntau = _mm_set1_pd(ntau);
#define d_APPLY_SMALL_DECL
#define d_APPLY_LARGE_DECL

#define d_DOT_FINISH_SMALL                      \
    vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp); \
    vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);

#define d_DOT_FINISH_LARGE                           \
    vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);     \
    vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1); \
    vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);       \
    vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);       \
    vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);      \
    vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);       \
    vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);     \
    vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

#define s_ELT real
#define s_V2 __m128
#define s_V4 __m256
#define s_W2 4
#define s_W4 8

#define s_THRES_SML 1.0842022e-19f
#define s_THRES_BIG 4.5035996e+15f
#define s_SCALE_SML 3.7778932e+22f
#define s_SCALE_BIG 1.3234890e-23f

#define s_LAMCH(s) slamch_(s)
#define s_SCAL aocl_blas_sscal

#define s_V2_SET1(x) _mm_set1_ps(x)
#define s_V4_SET1(x) _mm256_set1_ps(x)
#define s_V2_SETZERO() _mm_setzero_ps()
#define s_V4_SETZERO() _mm256_setzero_ps()
#define s_V2_LOADU(p) _mm_loadu_ps((const real *)(p))
#define s_V4_LOADU(p) _mm256_loadu_ps((const real *)(p))
#define s_V2_STOREU(p, v) _mm_storeu_ps((real *)(p), v)
#define s_V4_STOREU(p, v) _mm256_storeu_ps((real *)(p), v)
#define s_V2_LOAD_S(p) _mm_load_ss((const real *)(p))
#define s_V2_STORE_S(p, v) _mm_store_ss((real *)(p), v)
#define s_V2_LOAD_DUP(p) _mm_load_ss((const real *)(p))
#define s_V2_STORE_L(p, v) _mm_store_ss((real *)(p), v)
#define s_V2_LOAD_DUP_NORM(p) _mm_load1_ps((const real *)(p))
#define s_V2_MUL(a, b) _mm_mul_ps(a, b)
#define s_V4_MUL(a, b) _mm256_mul_ps(a, b)
#define s_V2_ADD(a, b) _mm_add_ps(a, b)
#define s_V4_ADD(a, b) _mm256_add_ps(a, b)
#define s_V2_ANDNOT(a, b) _mm_andnot_ps(a, b)
#define s_V4_ANDNOT(a, b) _mm256_andnot_ps(a, b)
#define s_V2_CMP(a, b, c) _mm_cmp_ps(a, b, c)
#define s_V4_CMP(a, b, c) _mm256_cmp_ps(a, b, c)
#define s_V2_OR(a, b) _mm_or_ps(a, b)
#define s_V4_OR(a, b) _mm256_or_ps(a, b)
#define s_V2_TESTZ(a, b) _mm_testz_ps(a, b)
#define s_V4_TESTZ(a, b) _mm256_testz_ps(a, b)
#define s_V2_BLENDV(a, b, m) _mm_blendv_ps(a, b, m)
#define s_V4_BLENDV(a, b, m) _mm256_blendv_ps(a, b, m)
#define s_V2_HSUM(v) ((v)[0] + (v)[1] + (v)[2] + (v)[3])
#define s_V4_HSUM(v) ((v)[0] + (v)[1] + (v)[2] + (v)[3] + (v)[4] + (v)[5] + (v)[6] + (v)[7])
#define s_TAIL_LOOP(x, lim) for(; (x) <= (lim); (x)++)

#define s_INIT_EXTRA
#define s_APPLY_NTAU_SETUP
#define s_APPLY_SMALL_DECL real sdot;
#define s_APPLY_LARGE_DECL \
    real sdot;             \
    __m128 hsum;

#define s_DOT_FINISH_SMALL                      \
    vd2_dtmp = _mm_hadd_ps(vd2_dtmp, vd2_dtmp); \
    vd2_dtmp = _mm_hadd_ps(vd2_dtmp, vd2_dtmp); \
    sdot = _mm_cvtss_f32(vd2_dtmp) * ntau;      \
    vd2_dtmp = _mm_set1_ps(sdot);

#define s_DOT_FINISH_LARGE                                   \
    hsum = _mm_add_ps(_mm256_castps256_ps128(vd4_dtmp),      \
                      _mm256_extractf128_ps(vd4_dtmp, 0x1)); \
    hsum = _mm_add_ps(hsum, vd2_dtmp);                       \
    hsum = _mm_hadd_ps(hsum, hsum);                          \
    hsum = _mm_hadd_ps(hsum, hsum);                          \
    sdot = _mm_cvtss_f32(hsum) * ntau;                       \
    vd4_dtmp = _mm256_set1_ps(sdot);                         \
    vd2_dtmp = _mm_set1_ps(sdot);

/* Declaration of local variables for QR Small */
#define FLA_GEQRF_INIT_GSMALL(P)                                                            \
    aocl_int64_t i, j, k;                                                                   \
    aocl_int64_t kcnt, slen;                                                                \
    aocl_int64_t acols, arows;                                                              \
    FLA_GEQRF_SMALL_CAT(P, _ELT) xnorm, vnorm, dtmp;                                        \
    FLA_GEQRF_SMALL_CAT(P, _ELT) fnorm, scale;                                              \
    FLA_GEQRF_SMALL_CAT(P, _ELT) med_sum, sml_sum, big_sum;                                 \
    FLA_GEQRF_SMALL_CAT(P, _ELT) alpha, beta, ntau;                                         \
    FLA_GEQRF_SMALL_CAT(P, _ELT) *v, *A, *ac;                                               \
                                                                                            \
    static TLS_CLASS_SPEC int r_once = 1;                                                   \
    static FLA_GEQRF_SMALL_CAT(P, _ELT) safmin, rsafmin;                                    \
    /* Constants chosen to minimize roundoff, */                                            \
    /* according to Blue's algorithm          */                                            \
    static FLA_GEQRF_SMALL_CAT(P, _ELT) thres_sml = FLA_GEQRF_SMALL_CAT(P, _THRES_SML);     \
    static FLA_GEQRF_SMALL_CAT(P, _ELT) thres_big = FLA_GEQRF_SMALL_CAT(P, _THRES_BIG);     \
    static FLA_GEQRF_SMALL_CAT(P, _ELT) scale_sml = FLA_GEQRF_SMALL_CAT(P, _SCALE_SML);     \
    static FLA_GEQRF_SMALL_CAT(P, _ELT) scale_big = FLA_GEQRF_SMALL_CAT(P, _SCALE_BIG);     \
                                                                                            \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_inp, vd2_abs_inp;                                       \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_sth, vd2_bth, vd2_sscl, vd2_bscl;                       \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_sinp, vd2_binp, vd2_minp;                               \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_smsk, vd2_bmsk, vd2_mmsk;                               \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_ssum, vd2_bsum, vd2_msum;                               \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_norm, vd2_vj1;                                          \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_dtmp, vd2_dtmp2, vd2_dtmp3;                             \
    FLA_GEQRF_SMALL_CAT(P, _INIT_EXTRA)                                                     \
    FLA_GEQRF_SMALL_CAT(P, _V2) vd2_zero = FLA_GEQRF_SMALL_CAT(P, _V2_SET1)(-0.0f);         \
                                                                                            \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_inp, vd4_abs_inp;                                       \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_sth, vd4_bth, vd4_sscl, vd4_bscl;                       \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_sinp, vd4_binp, vd4_minp;                               \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_smsk, vd4_bmsk, vd4_mmsk;                               \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_ssum, vd4_bsum, vd4_msum;                               \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_norm, vd4_vj;                                           \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_dtmp, vd4_dtmp2, vd4_dtmp3;                             \
    FLA_GEQRF_SMALL_CAT(P, _V4) vd4_zero = FLA_GEQRF_SMALL_CAT(P, _V4_SET1)(-0.0f);         \
                                                                                            \
    if(r_once)                                                                              \
    {                                                                                       \
        safmin = FLA_GEQRF_SMALL_CAT(P, _LAMCH)("S") / FLA_GEQRF_SMALL_CAT(P, _LAMCH)("E"); \
        rsafmin = 1. / safmin;                                                              \
        r_once = 0;                                                                         \
    }                                                                                       \
                                                                                            \
    vd2_sth = FLA_GEQRF_SMALL_CAT(P, _V2_SET1)(thres_sml);                                  \
    vd2_bth = FLA_GEQRF_SMALL_CAT(P, _V2_SET1)(thres_big);                                  \
    vd2_sscl = FLA_GEQRF_SMALL_CAT(P, _V2_SET1)(scale_sml);                                 \
    vd2_bscl = FLA_GEQRF_SMALL_CAT(P, _V2_SET1)(scale_big);

/* Combining three accumulators of norm to get final norm */
#define FLA_GEQRF_SMALL_GET_GNORM(P)                                     \
    fnorm = 0.0;                                                         \
    if(big_sum > 0.0)                                                    \
    {                                                                    \
        fnorm = big_sum;                                                 \
        if(med_sum > 0.0)                                                \
        {                                                                \
            fnorm = fnorm + (med_sum * scale_big) * scale_big;           \
        }                                                                \
        scale = scale_big;                                               \
    }                                                                    \
    else /* small sum must be non-zero */                                \
    {                                                                    \
        FLA_GEQRF_SMALL_CAT(P, _ELT) ymin, ymax;                         \
        if(med_sum > 0.0)                                                \
        {                                                                \
            med_sum = sqrt(med_sum);                                     \
            sml_sum = sqrt(sml_sum) / scale_sml;                         \
                                                                         \
            ymin = fla_min(med_sum, sml_sum);                            \
            ymax = fla_max(med_sum, sml_sum);                            \
                                                                         \
            scale = 1.0;                                                 \
            fnorm = ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax)); \
        }                                                                \
        else                                                             \
        {                                                                \
            scale = scale_sml;                                           \
            fnorm = sml_sum;                                             \
        }                                                                \
    }                                                                    \
    xnorm = sqrt(fnorm) / scale;

/* NORM computation using 256-bit AVX2 intrinsics */
#define FLA_GEQRF_SMALL_CALC_GNORM4(P, idx)                                             \
    /* load input and get its absolute values */                                        \
    vd4_inp = FLA_GEQRF_SMALL_CAT(P, _V4_LOADU)(&iptr[idx]);                            \
    vd4_abs_inp = FLA_GEQRF_SMALL_CAT(P, _V4_ANDNOT)(vd4_zero, vd4_inp);                \
                                                                                        \
    /* segregate input values into small, medium and big */                             \
    vd4_smsk = FLA_GEQRF_SMALL_CAT(P, _V4_CMP)(vd4_abs_inp, vd4_sth, _CMP_LT_OQ);       \
    vd4_bmsk = FLA_GEQRF_SMALL_CAT(P, _V4_CMP)(vd4_abs_inp, vd4_bth, _CMP_GT_OQ);       \
    vd4_mmsk = FLA_GEQRF_SMALL_CAT(P, _V4_OR)(vd4_smsk, vd4_bmsk);                      \
                                                                                        \
    /* if all inputs are in  medium range */                                            \
    if(FLA_GEQRF_SMALL_CAT(P, _V4_TESTZ)(vd4_mmsk, vd4_mmsk))                           \
    {                                                                                   \
        vd4_dtmp = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_inp, vd4_inp);                   \
        vd4_msum = FLA_GEQRF_SMALL_CAT(P, _V4_ADD)(vd4_msum, vd4_dtmp);                 \
    }                                                                                   \
    else /* for small and large inputs */                                               \
    {                                                                                   \
        has_outliers = 1;                                                               \
        vd4_sinp = FLA_GEQRF_SMALL_CAT(P, _V4_BLENDV)(vd4_zero, vd4_abs_inp, vd4_smsk); \
        vd4_binp = FLA_GEQRF_SMALL_CAT(P, _V4_BLENDV)(vd4_zero, vd4_abs_inp, vd4_bmsk); \
        vd4_minp = FLA_GEQRF_SMALL_CAT(P, _V4_BLENDV)(vd4_abs_inp, vd4_zero, vd4_mmsk); \
                                                                                        \
        /* scale, square and add as applicable */                                       \
        vd4_sinp = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_sinp, vd4_sscl);                 \
        vd4_binp = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_binp, vd4_bscl);                 \
                                                                                        \
        vd4_dtmp = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_minp, vd4_minp);                 \
        vd4_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_sinp, vd4_sinp);                \
        vd4_dtmp3 = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_binp, vd4_binp);                \
        vd4_msum = FLA_GEQRF_SMALL_CAT(P, _V4_ADD)(vd4_msum, vd4_dtmp);                 \
        vd4_ssum = FLA_GEQRF_SMALL_CAT(P, _V4_ADD)(vd4_ssum, vd4_dtmp2);                \
        vd4_bsum = FLA_GEQRF_SMALL_CAT(P, _V4_ADD)(vd4_bsum, vd4_dtmp3);                \
    }

/* NORM computation using 128-bit AVX intrinsics */
#define FLA_GEQRF_SMALL_CALC_GNORM2(P)                                                  \
    /* get absolute value of the vector input */                                        \
    vd2_abs_inp = FLA_GEQRF_SMALL_CAT(P, _V2_ANDNOT)(vd2_zero, vd2_inp);                \
                                                                                        \
    /* compute flags to detect out-of-range values */                                   \
    vd2_smsk = FLA_GEQRF_SMALL_CAT(P, _V2_CMP)(vd2_abs_inp, vd2_sth, _CMP_LT_OQ);       \
    vd2_bmsk = FLA_GEQRF_SMALL_CAT(P, _V2_CMP)(vd2_abs_inp, vd2_bth, _CMP_GT_OQ);       \
    vd2_mmsk = FLA_GEQRF_SMALL_CAT(P, _V2_OR)(vd2_smsk, vd2_bmsk);                      \
                                                                                        \
    /* if all inputs are in  medium range */                                            \
    if(FLA_GEQRF_SMALL_CAT(P, _V2_TESTZ)(vd2_mmsk, vd2_mmsk))                           \
    {                                                                                   \
        vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_inp, vd2_inp);                   \
        vd2_msum = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_msum, vd2_dtmp);                 \
    }                                                                                   \
    else /* for small and large inputs */                                               \
    {                                                                                   \
        has_outliers = 1;                                                               \
        vd2_sinp = FLA_GEQRF_SMALL_CAT(P, _V2_BLENDV)(vd2_zero, vd2_abs_inp, vd2_smsk); \
        vd2_binp = FLA_GEQRF_SMALL_CAT(P, _V2_BLENDV)(vd2_zero, vd2_abs_inp, vd2_bmsk); \
        vd2_minp = FLA_GEQRF_SMALL_CAT(P, _V2_BLENDV)(vd2_abs_inp, vd2_zero, vd2_mmsk); \
                                                                                        \
        /* scale, square and add as applicable */                                       \
        vd2_sinp = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_sinp, vd2_sscl);                 \
        vd2_binp = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_binp, vd2_bscl);                 \
                                                                                        \
        vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_sinp, vd2_sinp);                 \
        vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_binp, vd2_binp);                \
        vd2_dtmp3 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_minp, vd2_minp);                \
        vd2_ssum = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_ssum, vd2_dtmp);                 \
        vd2_bsum = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_bsum, vd2_dtmp2);                \
        vd2_msum = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_msum, vd2_dtmp3);                \
    }

#define FLA_LARF_GEN_GSMALL_COL(P, i, m, n, tau)                          \
    /* calculate norm of sub-diagonal elements in current column */       \
    med_sum = sml_sum = big_sum = 0.;                                     \
    vd2_msum = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                     \
    vd2_ssum = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                     \
    vd2_bsum = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                     \
                                                                          \
    for(j = 1; j <= (slen - (FLA_GEQRF_SMALL_CAT(P, _W2) - 1));           \
        j += FLA_GEQRF_SMALL_CAT(P, _W2))                                 \
    {                                                                     \
        vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&iptr[j]);            \
        FLA_GEQRF_SMALL_CALC_GNORM2(P);                                   \
    }                                                                     \
                                                                          \
    FLA_GEQRF_SMALL_CAT(P, _TAIL_LOOP)(j, slen)                           \
    {                                                                     \
        /* load input and get its absolute values */                      \
        vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&iptr[j]);           \
        FLA_GEQRF_SMALL_CALC_GNORM2(P);                                   \
    }                                                                     \
                                                                          \
    /* Get all the three sums */                                          \
    med_sum = FLA_GEQRF_SMALL_CAT(P, _V2_HSUM)(vd2_msum);                 \
    /* Combining outlier accumulators if non-zero */                      \
    if(has_outliers)                                                      \
    {                                                                     \
        sml_sum = FLA_GEQRF_SMALL_CAT(P, _V2_HSUM)(vd2_ssum);             \
        big_sum = FLA_GEQRF_SMALL_CAT(P, _V2_HSUM)(vd2_bsum);             \
        FLA_GEQRF_SMALL_GET_GNORM(P);                                     \
    }                                                                     \
    else                                                                  \
    {                                                                     \
        xnorm = sqrt(med_sum);                                            \
    }                                                                     \
                                                                          \
    /* Compute Householder Reflector parameters */                        \
    if(xnorm == 0.) /* Sub-diagonal elements are already zero */          \
    {                                                                     \
        tau[i] = 0.;                                                      \
        beta = 0.;                                                        \
    }                                                                     \
    else /* Non-zero sub-diagonal elements */                             \
    {                                                                     \
        /* Part 1: Compute Householder vector 'v' and tau */              \
                                                                          \
        v = iptr - 1;                                                     \
        alpha = v[1];                                                     \
        /* check for NAN */                                               \
        if(alpha != alpha || xnorm != xnorm)                              \
        {                                                                 \
            beta = alpha + xnorm;                                         \
        }                                                                 \
        else                                                              \
        {                                                                 \
            FLA_GEQRF_SMALL_CAT(P, _ELT) w, z;                            \
                                                                          \
            dtmp = f2c_abs(alpha);                                        \
            w = fla_max(dtmp, xnorm);                                     \
            z = fla_min(dtmp, xnorm);                                     \
                                                                          \
            z = z / w;                                                    \
            beta = w * sqrt(z * z + 1);                                   \
        }                                                                 \
        beta = (alpha >= 0.) ? -beta : beta;                              \
                                                                          \
        /* Scale-up the inputs for small norm */                          \
        for(kcnt = 0; (f2c_abs(beta) < safmin && kcnt <= 20); kcnt++)     \
        {                                                                 \
            FLA_GEQRF_SMALL_CAT(P, _SCAL)(&slen, &rsafmin, &v[2], &c__1); \
            beta = beta * rsafmin;                                        \
            alpha = alpha * rsafmin;                                      \
        }                                                                 \
                                                                          \
        /* Calculate tau and v */                                         \
        tau[i] = (beta - alpha) / beta;                                   \
        vnorm = 1. / (alpha - beta);                                      \
        /* Scale current column by norm to get v */                       \
        vd2_norm = FLA_GEQRF_SMALL_CAT(P, _V2_SET1)(vnorm);               \
                                                                          \
        /* Normalize using SIMD */                                        \
        for(j = 1; j <= (slen - (FLA_GEQRF_SMALL_CAT(P, _W2) - 1));       \
            j += FLA_GEQRF_SMALL_CAT(P, _W2))                             \
        {                                                                 \
            vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&v[j + 1]);       \
            vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_norm); \
            FLA_GEQRF_SMALL_CAT(P, _V2_STOREU)(&v[j + 1], vd2_vj1);       \
        }                                                                 \
        FLA_GEQRF_SMALL_CAT(P, _TAIL_LOOP)(j, slen)                       \
        {                                                                 \
            vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_DUP)(&v[j + 1]);    \
            vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_norm); \
            FLA_GEQRF_SMALL_CAT(P, _V2_STORE_L)(&v[j + 1], vd2_vj1);      \
        }                                                                 \
        /* Scale-down beta */                                             \
        for(; kcnt >= 1; kcnt--)                                          \
        {                                                                 \
            beta = beta * safmin;                                         \
        }                                                                 \
    }

#define FLA_LARF_APPLY_GSMALL_COL(P, i, m, n, r, ldr, tau)                       \
    if(xnorm != 0.) /* Sub-diagonal elements are already zero */                 \
    {                                                                            \
        /* Part 2: Apply the Householder rotation              */                \
        /* on the rest of the matrix                           */                \
        /*    A = A - tau * v * v**T * A                       */                \
        /*      = A - v * tau * (A**T * v)**T                  */                \
                                                                                 \
        FLA_GEQRF_SMALL_CAT(P, _APPLY_SMALL_DECL)                                \
        A = &r[i + (i + 1) * *ldr];                                              \
        arows = *m - i + 1;                                                      \
        acols = *n - i;                                                          \
        v[1] = 1.;                                                               \
        ntau = -tau[i];                                                          \
        FLA_GEQRF_SMALL_CAT(P, _APPLY_NTAU_SETUP)                                \
                                                                                 \
        /* Compute A**T * v */                                                   \
        for(j = 1; j <= acols; j++) /* for every column c_A of A */              \
        {                                                                        \
            ac = &A[(j - 1) * *ldr - 1];                                         \
                                                                                 \
            /* Compute tmp = c_A**T . v */                                       \
            vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                    \
            for(k = 1; k <= (arows - (FLA_GEQRF_SMALL_CAT(P, _W2) - 1));         \
                k += FLA_GEQRF_SMALL_CAT(P, _W2))                                \
            {                                                                    \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&ac[k]);             \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&v[k]);              \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_inp, vd2_vj1);   \
                vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_dtmp, vd2_dtmp2); \
            }                                                                    \
            FLA_GEQRF_SMALL_CAT(P, _TAIL_LOOP)(k, arows)                         \
            {                                                                    \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&ac[k]);            \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&v[k]);             \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_inp, vd2_vj1);   \
                vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_dtmp, vd2_dtmp2); \
            }                                                                    \
                                                                                 \
            /* Compute tmp = -tau * (c_A**T . v) and broadcast */                \
            FLA_GEQRF_SMALL_CAT(P, _DOT_FINISH_SMALL)                            \
                                                                                 \
            /* Compute c_A + tmp * v */                                          \
            for(k = 1; k <= (arows - (FLA_GEQRF_SMALL_CAT(P, _W2) - 1));         \
                k += FLA_GEQRF_SMALL_CAT(P, _W2))                                \
            {                                                                    \
                /* load column elements of c_A and v */                          \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&ac[k]);             \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&v[k]);              \
                                                                                 \
                /* mul by dtmp, add and store */                                 \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_dtmp);  \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_inp, vd2_dtmp2);   \
                FLA_GEQRF_SMALL_CAT(P, _V2_STOREU)(&ac[k], vd2_inp);             \
            }                                                                    \
            FLA_GEQRF_SMALL_CAT(P, _TAIL_LOOP)(k, arows)                         \
            {                                                                    \
                /* load single remaining element from c_A and v */               \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&ac[k]);            \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&v[k]);             \
                                                                                 \
                /* multiply with tau and store */                                \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_dtmp);  \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_inp, vd2_dtmp2);   \
                FLA_GEQRF_SMALL_CAT(P, _V2_STORE_S)(&ac[k], vd2_inp);            \
            }                                                                    \
        }                                                                        \
        v[1] = beta;                                                             \
    }

#define FLA_LARF_GEN_GLARGE_COL(P, i, m, n, tau)                              \
    /* calculate norm of sub-diagonal elements in current column */           \
    med_sum = sml_sum = big_sum = 0.;                                         \
    vd4_sth = FLA_GEQRF_SMALL_CAT(P, _V4_SET1)(thres_sml);                    \
    vd4_bth = FLA_GEQRF_SMALL_CAT(P, _V4_SET1)(thres_big);                    \
    vd4_sscl = FLA_GEQRF_SMALL_CAT(P, _V4_SET1)(scale_sml);                   \
    vd4_bscl = FLA_GEQRF_SMALL_CAT(P, _V4_SET1)(scale_big);                   \
                                                                              \
    vd4_msum = FLA_GEQRF_SMALL_CAT(P, _V4_SETZERO)();                         \
    vd4_ssum = FLA_GEQRF_SMALL_CAT(P, _V4_SETZERO)();                         \
    vd4_bsum = FLA_GEQRF_SMALL_CAT(P, _V4_SETZERO)();                         \
                                                                              \
    for(j = 1; j <= (slen - (FLA_GEQRF_SMALL_CAT(P, _W4) - 1));               \
        j += FLA_GEQRF_SMALL_CAT(P, _W4))                                     \
    {                                                                         \
        FLA_GEQRF_SMALL_CALC_GNORM4(P, j);                                    \
    }                                                                         \
                                                                              \
    if(j <= slen)                                                             \
    { /* process remaining iterations */                                      \
        vd2_msum = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                     \
        vd2_ssum = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                     \
        vd2_bsum = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                     \
                                                                              \
        /* process remaining inputs per iteration */                          \
        for(; j <= slen; j++)                                                 \
        {                                                                     \
            vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_DUP_NORM)(&iptr[j]);    \
            FLA_GEQRF_SMALL_CALC_GNORM2(P);                                   \
        }                                                                     \
        /* Get all the three sums */                                          \
        med_sum = FLA_GEQRF_SMALL_CAT(P, _V4_HSUM)(vd4_msum) + vd2_msum[0];   \
        sml_sum = FLA_GEQRF_SMALL_CAT(P, _V4_HSUM)(vd4_ssum) + vd2_ssum[0];   \
        big_sum = FLA_GEQRF_SMALL_CAT(P, _V4_HSUM)(vd4_bsum) + vd2_bsum[0];   \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        /* Get all the three sums in case of no remaining iterations */       \
        med_sum = FLA_GEQRF_SMALL_CAT(P, _V4_HSUM)(vd4_msum);                 \
        sml_sum = FLA_GEQRF_SMALL_CAT(P, _V4_HSUM)(vd4_ssum);                 \
        big_sum = FLA_GEQRF_SMALL_CAT(P, _V4_HSUM)(vd4_bsum);                 \
    }                                                                         \
                                                                              \
    /* Combining outlier accumulators if non-zero */                          \
    if(has_outliers)                                                          \
    {                                                                         \
        FLA_GEQRF_SMALL_GET_GNORM(P);                                         \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        xnorm = sqrt(med_sum);                                                \
    }                                                                         \
                                                                              \
    /* Compute Householder Reflector parameters */                            \
    if(xnorm == 0.) /* Sub-diagonal elements are already zero */              \
    {                                                                         \
        tau[i] = 0.;                                                          \
        beta = 0.;                                                            \
    }                                                                         \
    else /* Non-zero sub-diagonal elements */                                 \
    {                                                                         \
        /* Part 1: Compute Householder vector 'v' and tau */                  \
                                                                              \
        v = iptr - 1;                                                         \
        alpha = v[1];                                                         \
                                                                              \
        /* Compute Householder rotated vector */                              \
        if(alpha != alpha || xnorm != xnorm) /* check for NAN */              \
        {                                                                     \
            beta = alpha + xnorm;                                             \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            FLA_GEQRF_SMALL_CAT(P, _ELT) w, z;                                \
                                                                              \
            dtmp = f2c_abs(alpha);                                            \
            w = fla_max(dtmp, xnorm);                                         \
            z = fla_min(dtmp, xnorm);                                         \
                                                                              \
            z = z / w;                                                        \
            beta = w * sqrt(z * z + 1);                                       \
        }                                                                     \
        beta = (alpha >= 0.) ? -beta : beta;                                  \
                                                                              \
        /* Scale-up the inputs for small norm */                              \
        for(kcnt = 0; (f2c_abs(beta) < safmin && kcnt <= 20); kcnt++)         \
        {                                                                     \
            FLA_GEQRF_SMALL_CAT(P, _SCAL)(&slen, &rsafmin, &v[2], &c__1);     \
            beta = beta * rsafmin;                                            \
            alpha = alpha * rsafmin;                                          \
        }                                                                     \
                                                                              \
        /* Calculate tau and v */                                             \
        tau[i] = (beta - alpha) / beta;                                       \
        vnorm = 1. / (alpha - beta);                                          \
        /* Scale current column by norm to get v */                           \
        vd4_norm = FLA_GEQRF_SMALL_CAT(P, _V4_SET1)(vnorm);                   \
                                                                              \
        /* Normalize using SIMD */                                            \
        for(j = 1; j <= (slen - (FLA_GEQRF_SMALL_CAT(P, _W4) - 1));           \
            j += FLA_GEQRF_SMALL_CAT(P, _W4))                                 \
        {                                                                     \
            vd4_vj = FLA_GEQRF_SMALL_CAT(P, _V4_LOADU)(&v[j + 1]);            \
            vd4_vj = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_vj, vd4_norm);       \
            FLA_GEQRF_SMALL_CAT(P, _V4_STOREU)(&v[j + 1], vd4_vj);            \
        }                                                                     \
        /* Remaining iterations through 128-bit SIMD */                       \
        if(j <= slen)                                                         \
        {                                                                     \
            vd2_norm = FLA_GEQRF_SMALL_CAT(P, _V2_SET1)(vnorm);               \
            for(; j <= (slen - (FLA_GEQRF_SMALL_CAT(P, _W2) - 1));            \
                j += FLA_GEQRF_SMALL_CAT(P, _W2))                             \
            {                                                                 \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&v[j + 1]);       \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_norm); \
                FLA_GEQRF_SMALL_CAT(P, _V2_STOREU)(&v[j + 1], vd2_vj1);       \
            }                                                                 \
            FLA_GEQRF_SMALL_CAT(P, _TAIL_LOOP)(j, slen)                       \
            {                                                                 \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_DUP)(&v[j + 1]);    \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_norm); \
                FLA_GEQRF_SMALL_CAT(P, _V2_STORE_L)(&v[j + 1], vd2_vj1);      \
            }                                                                 \
        }                                                                     \
        /* Scale-down beta */                                                 \
        for(; kcnt >= 1; kcnt--)                                              \
        {                                                                     \
            beta = beta * safmin;                                             \
        }                                                                     \
    }

#define FLA_LARF_APPLY_GLARGE_COL(P, i, m, n, r, ldr, tau)                       \
    if(xnorm != 0.) /* Sub-diagonal elements are already zero */                 \
    {                                                                            \
        /* Part 2: Apply the Householder rotation              */                \
        /* on the rest of the matrix                           */                \
        /*    A = A - tau * v * v**T * A                       */                \
        /*      = A - v * tau * (A**T * v)**T                  */                \
                                                                                 \
        FLA_GEQRF_SMALL_CAT(P, _APPLY_LARGE_DECL)                                \
        A = &r[i + (i + 1) * *ldr];                                              \
        arows = *m - i + 1;                                                      \
        acols = *n - i;                                                          \
        v[1] = 1.;                                                               \
        ntau = -tau[i];                                                          \
        FLA_GEQRF_SMALL_CAT(P, _APPLY_NTAU_SETUP)                                \
                                                                                 \
        /* Compute A**T * v */                                                   \
        for(j = 1; j <= acols; j++) /* for every column c_A of A */              \
        {                                                                        \
            ac = &A[(j - 1) * *ldr - 1];                                         \
            vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_SETZERO)();                    \
            vd4_dtmp = FLA_GEQRF_SMALL_CAT(P, _V4_SETZERO)();                    \
                                                                                 \
            /* Compute tmp = c_A**T . v */                                       \
            for(k = 1; k <= (arows - (FLA_GEQRF_SMALL_CAT(P, _W4) - 1));         \
                k += FLA_GEQRF_SMALL_CAT(P, _W4))                                \
            {                                                                    \
                /* load column elements of A and v */                            \
                vd4_inp = FLA_GEQRF_SMALL_CAT(P, _V4_LOADU)(&ac[k]);             \
                vd4_vj = FLA_GEQRF_SMALL_CAT(P, _V4_LOADU)(&v[k]);               \
                                                                                 \
                /* take dot product */                                           \
                vd4_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_inp, vd4_vj);    \
                vd4_dtmp = FLA_GEQRF_SMALL_CAT(P, _V4_ADD)(vd4_dtmp, vd4_dtmp2); \
            }                                                                    \
            if(k <= (arows - (FLA_GEQRF_SMALL_CAT(P, _W2) - 1)))                 \
            {                                                                    \
                /* load column elements of A and v */                            \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&ac[k]);             \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&v[k]);              \
                                                                                 \
                /* take dot product */                                           \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_inp, vd2_vj1);   \
                vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_dtmp, vd2_dtmp2); \
                k += FLA_GEQRF_SMALL_CAT(P, _W2);                                \
            }                                                                    \
            FLA_GEQRF_SMALL_CAT(P, _TAIL_LOOP)(k, arows)                         \
            {                                                                    \
                /* load single remaining element from c_A and v */               \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&ac[k]);            \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&v[k]);             \
                                                                                 \
                /* take dot product */                                           \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_inp, vd2_vj1);   \
                vd2_dtmp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_dtmp, vd2_dtmp2); \
            }                                                                    \
            /* Horizontal add of dtmp */                                         \
            FLA_GEQRF_SMALL_CAT(P, _DOT_FINISH_LARGE)                            \
                                                                                 \
            /* Compute c_A + tmp * v */                                          \
            for(k = 1; k <= (arows - (FLA_GEQRF_SMALL_CAT(P, _W4) - 1));         \
                k += FLA_GEQRF_SMALL_CAT(P, _W4))                                \
            {                                                                    \
                /* load column elements of c_A and v */                          \
                vd4_inp = FLA_GEQRF_SMALL_CAT(P, _V4_LOADU)(&ac[k]);             \
                vd4_vj = FLA_GEQRF_SMALL_CAT(P, _V4_LOADU)(&v[k]);               \
                                                                                 \
                /* mul by dtmp, add and store */                                 \
                vd4_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V4_MUL)(vd4_vj, vd4_dtmp);   \
                vd4_inp = FLA_GEQRF_SMALL_CAT(P, _V4_ADD)(vd4_inp, vd4_dtmp2);   \
                FLA_GEQRF_SMALL_CAT(P, _V4_STOREU)(&ac[k], vd4_inp);             \
            }                                                                    \
            if(k <= (arows - (FLA_GEQRF_SMALL_CAT(P, _W2) - 1)))                 \
            {                                                                    \
                /* load column elements of c_A and v */                          \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&ac[k]);             \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOADU)(&v[k]);              \
                                                                                 \
                /* mul by dtmp, add and store */                                 \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_dtmp);  \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_inp, vd2_dtmp2);   \
                FLA_GEQRF_SMALL_CAT(P, _V2_STOREU)(&ac[k], vd2_inp);             \
                k += FLA_GEQRF_SMALL_CAT(P, _W2);                                \
            }                                                                    \
            FLA_GEQRF_SMALL_CAT(P, _TAIL_LOOP)(k, arows)                         \
            {                                                                    \
                /* load single remaining element from c_A and v */               \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&ac[k]);            \
                vd2_vj1 = FLA_GEQRF_SMALL_CAT(P, _V2_LOAD_S)(&v[k]);             \
                                                                                 \
                /* mul by dtmp, add and store */                                 \
                vd2_dtmp2 = FLA_GEQRF_SMALL_CAT(P, _V2_MUL)(vd2_vj1, vd2_dtmp);  \
                vd2_inp = FLA_GEQRF_SMALL_CAT(P, _V2_ADD)(vd2_inp, vd2_dtmp2);   \
                FLA_GEQRF_SMALL_CAT(P, _V2_STORE_L)(&ac[k], vd2_inp);            \
            }                                                                    \
        }                                                                        \
        v[1] = beta;                                                             \
    }

/* Double precision */
#define FLA_GEQRF_INIT_DSMALL() FLA_GEQRF_INIT_GSMALL(d)
#define FLA_LARF_GEN_DSMALL_COL(i, m, n, tau) FLA_LARF_GEN_GSMALL_COL(d, i, m, n, tau)
#define FLA_LARF_APPLY_DSMALL_COL(i, m, n, r, ldr, tau) \
    FLA_LARF_APPLY_GSMALL_COL(d, i, m, n, r, ldr, tau)
#define FLA_LARF_GEN_DLARGE_COL(i, m, n, tau) FLA_LARF_GEN_GLARGE_COL(d, i, m, n, tau)
#define FLA_LARF_APPLY_DLARGE_COL(i, m, n, r, ldr, tau) \
    FLA_LARF_APPLY_GLARGE_COL(d, i, m, n, r, ldr, tau)

/* Single precision */
#define FLA_GEQRF_SLARGE_MIN_LEN s_W4
#define FLA_GEQRF_INIT_SSMALL() FLA_GEQRF_INIT_GSMALL(s)
#define FLA_LARF_GEN_SSMALL_COL(i, m, n, tau) FLA_LARF_GEN_GSMALL_COL(s, i, m, n, tau)
#define FLA_LARF_APPLY_SSMALL_COL(i, m, n, r, ldr, tau) \
    FLA_LARF_APPLY_GSMALL_COL(s, i, m, n, r, ldr, tau)
#define FLA_LARF_GEN_SLARGE_COL(i, m, n, tau) FLA_LARF_GEN_GLARGE_COL(s, i, m, n, tau)
#define FLA_LARF_APPLY_SLARGE_COL(i, m, n, r, ldr, tau) \
    FLA_LARF_APPLY_GLARGE_COL(s, i, m, n, r, ldr, tau)

#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_GEQRF_SMALL_AVX2_KERNEL_H */
