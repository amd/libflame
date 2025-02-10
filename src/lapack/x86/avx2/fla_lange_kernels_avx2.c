/******************************************************************************
 * * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

#define invoke_mm256_loadu_px(x, ...) _mm256_loadu_p##x(__VA_ARGS__)
#define invoke_mm256_permute2f128_px(x, ...) _mm256_permute2f128_p##x(__VA_ARGS__)
#define invoke_mm256_permute_px(x, ...) _mm256_permute_p##x(__VA_ARGS__)
#define invoke_mm256_max_px(x, ...) _mm256_max_p##x(__VA_ARGS__)
#define invoke_mm256_andnot_px(x, ...) _mm256_andnot_p##x(__VA_ARGS__)
#define invoke_mm256_set1_px(x, ...) _mm256_set1_p##x(__VA_ARGS__)
#define invoke_mm256_shuffle_px(x, ...) _mm256_shuffle_p##x(__VA_ARGS__)
#define invoke_mm256_mul_px(x, ...) _mm256_mul_p##x(__VA_ARGS__)
#define invoke_mm256_add_px(x, ...) _mm256_add_p##x(__VA_ARGS__)
#define invoke_mm256_castpx256_px128(x, ...) _mm256_castp##x##256_p##x##128(__VA_ARGS__)
#define invoke_mm256_extractf128_px(x, ...) _mm256_extractf128_p##x(__VA_ARGS__)
#define invoke_mm256_cmp_px(x, ...) _mm256_cmp_p##x(__VA_ARGS__)
#define invoke_mm256_or_px(x, ...) _mm256_or_p##x(__VA_ARGS__)
#define invoke_mm256_movemask_px(x, ...) _mm256_movemask_p##x(__VA_ARGS__)

#define invoke_mm_loadu_px(x, ...) _mm_loadu_p##x(__VA_ARGS__)
#define invoke_mm_max_px(x, ...) _mm_max_p##x(__VA_ARGS__)
#define invoke_mm_andnot_px(x, ...) _mm_andnot_p##x(__VA_ARGS__)
#define invoke_mm_set1_px(x, ...) _mm_set1_p##x(__VA_ARGS__)
#define invoke_mm_shuffle_px(x, ...) _mm_shuffle_p##x(__VA_ARGS__)
#define invoke_mm_mul_px(x, ...) _mm_mul_p##x(__VA_ARGS__)
#define invoke_mm_add_px(x, ...) _mm_add_p##x(__VA_ARGS__)
#define invoke_mm_permute_px(x, ...) _mm_permute_p##x(__VA_ARGS__)
#define invoke_mm_cmp_px(x, ...) _mm_cmp_p##x(__VA_ARGS__)
#define invoke_mm_or_px(x, ...) _mm_or_p##x(__VA_ARGS__)
#define invoke_mm_movemask_px(x, ...) _mm_movemask_p##x(__VA_ARGS__)

#define invoke_sqrts(data) sqrtf(data)
#define invoke_sqrtd(data) sqrt(data)
#define invoke_sqrtx(x, data) invoke_sqrt##x(data)

#define invoke_mm_cvtss(...) _mm_cvtss_f32(__VA_ARGS__)
#define invoke_mm_cvtsd(...) _mm_cvtsd_f64(__VA_ARGS__)
#define invoke_mm_cvtsx(x, ...) invoke_mm_cvts##x(__VA_ARGS__)

#define type_m256s __m256
#define type_m256d __m256d
#define type_mm256x(x) type_m256##x

#define type_m128s __m128
#define type_m128d __m128d
#define type_mm128x(x) type_m128##x

#define fla_get_max_sabs_element_real_avx2_max_reduce        \
    vec1_128 = invoke_mm256_castpx256_px128(s, vmax);        \
    vec2_128 = invoke_mm256_extractf128_px(s, vmax, 1);      \
    vec1_128 = invoke_mm_max_px(s, vec1_128, vec2_128);      \
    vmax128 = invoke_mm_max_px(s, vec1_128, vmax128);        \
    vec1_128 = invoke_mm_permute_px(s, vmax128, 0b01001110); \
    vmax128 = invoke_mm_max_px(s, vec1_128, vmax128);        \
    vec2_128 = invoke_mm_permute_px(s, vmax128, 0b10110001); \
    vmax128 = invoke_mm_max_px(s, vec2_128, vmax128);

#define fla_get_max_dabs_element_real_avx2_max_reduce   \
    vec1_128 = invoke_mm256_castpx256_px128(d, vmax);   \
    vec2_128 = invoke_mm256_extractf128_px(d, vmax, 1); \
    vec1_128 = invoke_mm_max_px(d, vec1_128, vec2_128); \
    vmax128 = invoke_mm_max_px(d, vec1_128, vmax128);   \
    vec1_128 = invoke_mm_permute_px(d, vmax128, 0b01);  \
    vmax128 = invoke_mm_max_px(d, vec1_128, vmax128);

#define fla_get_max_xabs_element_real_avx2_max_reduce(x) \
    fla_get_max_##x##abs_element_real_avx2_max_reduce

/* x is the first letter for realtype: s or d */
#define fla_get_max_xabs_element_vector_real_avx2(x, realtype, m, a, a_dim, nx, stepsize, n128)  \
    type_mm256x(x) vec1, vec2, vec3;                                                             \
    type_mm128x(x) vec1_128, vec2_128;                                                           \
    type_mm256x(x) neg_zero = invoke_mm256_set1_px(x, -0.0f);                                    \
    type_mm128x(x) neg_zero128 = invoke_mm_set1_px(x, -0.0f);                                    \
    type_mm256x(x) vmax = {0};                                                                   \
    type_mm128x(x) vmax128 = {0};                                                                \
    type_mm256x(x) nan_flag = invoke_mm256_set1_px(x, 0);                                        \
    type_mm128x(x) nan_flag128 = invoke_mm_set1_px(x, 0);                                        \
    integer i__, lasti, found_nan;                                                               \
    realtype max_value, temp;                                                                    \
    max_value = 0.0f;                                                                            \
    lasti = m - (m % stepsize);                                                                  \
    /* Get maximum upto size lasti */                                                            \
    for(i__ = 1; i__ <= lasti; i__ += stepsize)                                                  \
    {                                                                                            \
        vec1 = invoke_mm256_loadu_px(x, &a[i__ + a_dim]);                                        \
        vec2 = invoke_mm256_loadu_px(x, &a[i__ + nx + a_dim]);                                   \
        /* Check for NaN */                                                                      \
        nan_flag                                                                                 \
            = invoke_mm256_or_px(x, nan_flag, invoke_mm256_cmp_px(x, vec1, vec2, _CMP_UNORD_Q)); \
        vec1 = invoke_mm256_andnot_px(x, neg_zero, vec1);                                        \
        vec2 = invoke_mm256_andnot_px(x, neg_zero, vec2);                                        \
        vec3 = invoke_mm256_max_px(x, vec1, vec2);                                               \
        vmax = invoke_mm256_max_px(x, vec3, vmax);                                               \
    }                                                                                            \
    if((m - lasti) >= nx)                                                                        \
    {                                                                                            \
        vec1 = invoke_mm256_loadu_px(x, &a[lasti + 1 + a_dim]);                                  \
        nan_flag                                                                                 \
            = invoke_mm256_or_px(x, nan_flag, invoke_mm256_cmp_px(x, vec1, vec1, _CMP_UNORD_Q)); \
        vec1 = invoke_mm256_andnot_px(x, neg_zero, vec1);                                        \
        vmax = invoke_mm256_max_px(x, vec1, vmax);                                               \
        lasti += nx;                                                                             \
    }                                                                                            \
    if((m - lasti) >= n128)                                                                      \
    {                                                                                            \
        vec1_128 = invoke_mm_loadu_px(x, &a[lasti + 1 + a_dim]);                                 \
        nan_flag128 = invoke_mm_or_px(x, nan_flag128,                                            \
                                      invoke_mm_cmp_px(x, vec1_128, vec1_128, _CMP_UNORD_Q));    \
        vec1_128 = invoke_mm_andnot_px(x, neg_zero128, vec1_128);                                \
        vmax128 = invoke_mm_max_px(x, vec1_128, vmax128);                                        \
        lasti += n128;                                                                           \
    }                                                                                            \
    found_nan = invoke_mm256_movemask_px(x, nan_flag) | invoke_mm_movemask_px(x, nan_flag128);   \
    if(found_nan)                                                                                \
    {                                                                                            \
        max_value = NAN;                                                                         \
    }                                                                                            \
    else                                                                                         \
    {                                                                                            \
        fla_get_max_xabs_element_real_avx2_max_reduce(x) /* Maximum value */                     \
            max_value                                                                            \
            = invoke_mm_cvtsx(x, vmax128);                                                       \
        for(i__ = lasti + 1; i__ <= m; i__++)                                                    \
        {                                                                                        \
            temp = f2c_abs(a[i__ + a_dim]);                                                      \
            if(max_value < temp || temp != temp)                                                 \
            {                                                                                    \
                max_value = temp;                                                                \
            }                                                                                    \
        }                                                                                        \
    }

#define fla_get_max_sabs_element_complex_mm256_shuffle_flag_real 0b10001000
#define fla_get_max_sabs_element_complex_mm256_shuffle_flag_imag 0b11011101

#define fla_get_max_dabs_element_complex_mm256_shuffle_flag_real 0b0000
#define fla_get_max_dabs_element_complex_mm256_shuffle_flag_imag 0b1111

#define fla_get_max_sabs_element_complex_mm256_shuffle_ri_flag 0b10110001
#define fla_get_max_dabs_element_complex_mm256_shuffle_ri_flag 0b0101

#define fla_get_max_sabs_element_complex_mm128_shuffle_ri_flag 0b10110001
#define fla_get_max_dabs_element_complex_mm128_shuffle_ri_flag 0b01

#define fla_get_max_xabs_element_complex_mm256_shuffle_flag_real(x) \
    fla_get_max_##x##abs_element_complex_mm256_shuffle_flag_real

#define fla_get_max_xabs_element_complex_mm256_shuffle_flag_imag(x) \
    fla_get_max_##x##abs_element_complex_mm256_shuffle_flag_imag

#define fla_get_max_xabs_element_complex_mm256_shuffle_ri_flag(x) \
    fla_get_max_##x##abs_element_complex_mm256_shuffle_ri_flag

#define fla_get_max_xabs_element_complex_mm128_shuffle_ri_flag(x) \
    fla_get_max_##x##abs_element_complex_mm128_shuffle_ri_flag

#define fla_get_max_sabs_element_vector_complex_avx2_load256(target_vec, idx) \
    target_vec##_128s = _mm_loadu_ps((real *)&a[(idx) + a_dim]);              \
    target_vec = _mm256_cvtps_pd(target_vec##_128s);

#define fla_get_max_dabs_element_vector_complex_avx2_load256(target_vec, idx) \
    target_vec = _mm256_loadu_pd((doublereal *)&a[(idx) + a_dim]);

#define fla_get_max_xabs_element_vector_complex_avx2_load256(x, target_vec, idx) \
    fla_get_max_##x##abs_element_vector_complex_avx2_load256(target_vec, idx)

#define fla_get_max_sabs_element_vector_complex_avx2_load128(target_vec, idx)          \
    target_vec##s = _mm_castsi128_ps(_mm_loadu_si64((doublereal *)&a[(idx) + a_dim])); \
    target_vec = _mm_cvtps_pd(target_vec##s);

#define fla_get_max_dabs_element_vector_complex_avx2_load128(target_vec, idx) \
    target_vec = _mm_loadu_pd((doublereal *)&a[(idx) + a_dim]);

#define fla_get_max_xabs_element_vector_complex_avx2_load128(x, target_vec, idx) \
    fla_get_max_##x##abs_element_vector_complex_avx2_load128(target_vec, idx)

#define fla_get_max_sabs_element_vector_complex_avx2_define_extra_vars \
    __m128 vec1_128s, vec2_128s, vec5_128s, vec6_128s;

#define fla_get_max_dabs_element_vector_complex_avx2_define_extra_vars

#define fla_get_max_xabs_element_vector_complex_avx2_define_extra_vars(x) \
    fla_get_max_##x##abs_element_vector_complex_avx2_define_extra_vars

#define fla_get_max_xabs_element_vector_complex_avx2(x, m, a, a_dim)                           \
    __m256d vec1, vec2, vec3, vec4;                                                            \
    __m256d vec5, vec6, vec7, vec8;                                                            \
    __m128d vec1_128, vec2_128;                                                                \
    fla_get_max_xabs_element_vector_complex_avx2_define_extra_vars(x);                         \
    __m256d vmax = {0};                                                                        \
    __m128d vmax128 = {0};                                                                     \
    __m256d nan_flag = _mm256_set1_pd(0.0);                                                    \
    __m128d nan_flag128 = _mm_set1_pd(0.0);                                                    \
    doublereal max_value;                                                                      \
    integer i__, lasti, found_nan;                                                             \
    max_value = 0.0f;                                                                          \
    lasti = m - (m % 8);                                                                       \
    /* Get maximum upto size lasti */                                                          \
    for(i__ = 1; i__ <= lasti; i__ += 8)                                                       \
    {                                                                                          \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec1, i__);                    \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec2, i__ + 2);                \
        /* Check for NaN */                                                                    \
        nan_flag = _mm256_or_pd(nan_flag, _mm256_cmp_pd(vec1, vec2, _CMP_UNORD_Q));            \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                               \
        vec3 = _mm256_shuffle_pd(vec1, vec2, 0b0000);                                          \
        vec4 = _mm256_shuffle_pd(vec1, vec2, 0b1111);                                          \
        /* calculating r^2 and i^r */                                                          \
        vec3 = _mm256_mul_pd(vec3, vec3);                                                      \
        vec4 = _mm256_mul_pd(vec4, vec4);                                                      \
        /* vec3 = vec3 + vec4 */                                                               \
        vec3 = _mm256_add_pd(vec3, vec4);                                                      \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec5, i__ + 4);                \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec6, i__ + 6);                \
        /* Check for NaN */                                                                    \
        nan_flag = _mm256_or_pd(nan_flag, _mm256_cmp_pd(vec5, vec6, _CMP_UNORD_Q));            \
        /* vec7 = r = real(a[nx:2*nx]); vec8 = i = image(a[nx:2*nx])*/                         \
        vec7 = _mm256_shuffle_pd(vec5, vec6, 0b0000);                                          \
        vec8 = _mm256_shuffle_pd(vec5, vec6, 0b1111);                                          \
        /* calculating r^2 and i^r */                                                          \
        vec7 = _mm256_mul_pd(vec7, vec7);                                                      \
        vec8 = _mm256_mul_pd(vec8, vec8);                                                      \
        /* vec7 = vec7 + vec8 */                                                               \
        vec7 = _mm256_add_pd(vec7, vec8);                                                      \
        /* vec1 = max(vec3, vec7) */                                                           \
        vec1 = _mm256_max_pd(vec3, vec7);                                                      \
        vmax = _mm256_max_pd(vec1, vmax);                                                      \
    }                                                                                          \
    if((m - lasti) >= 4)                                                                       \
    {                                                                                          \
        /* Permute and find the maximum element*/                                              \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec1, lasti + 1);              \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec2, lasti + 3);              \
        /* Check for NaN */                                                                    \
        nan_flag = _mm256_or_pd(nan_flag, _mm256_cmp_pd(vec1, vec2, _CMP_UNORD_Q));            \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                               \
        vec3 = _mm256_shuffle_pd(vec1, vec2, 0b0000);                                          \
        vec4 = _mm256_shuffle_pd(vec1, vec2, 0b1111);                                          \
        /* calculating r^2 and i^r */                                                          \
        vec3 = _mm256_mul_pd(vec3, vec3);                                                      \
        vec4 = _mm256_mul_pd(vec4, vec4);                                                      \
        /* vec3 = vec3 + vec4 */                                                               \
        vec3 = _mm256_add_pd(vec3, vec4);                                                      \
        vmax = _mm256_max_pd(vec3, vmax);                                                      \
        lasti += 4;                                                                            \
    }                                                                                          \
    if((m - lasti) >= 2)                                                                       \
    {                                                                                          \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec1, lasti + 1);              \
        /* Check for NaN */                                                                    \
        nan_flag = _mm256_or_pd(nan_flag, _mm256_cmp_pd(vec1, vec1, _CMP_UNORD_Q));            \
        vec1 = _mm256_mul_pd(vec1, vec1);                                                      \
        vec2 = _mm256_permute_pd(vec1, 0b0101);                                                \
        vec1 = _mm256_add_pd(vec1, vec2);                                                      \
        vmax = _mm256_max_pd(vec1, vmax);                                                      \
        lasti += 2;                                                                            \
    }                                                                                          \
    if((m - lasti) == 1)                                                                       \
    {                                                                                          \
        fla_get_max_xabs_element_vector_complex_avx2_load128(x, vec1_128, m);                  \
        /* Check for NaN */                                                                    \
        nan_flag128 = _mm_or_pd(nan_flag128, _mm_cmp_pd(vec1_128, vec1_128, _CMP_UNORD_Q));    \
        vec1_128 = _mm_mul_pd(vec1_128, vec1_128);                                             \
        vec2_128 = _mm_permute_pd(vec1_128, 0b01);                                             \
        vmax128 = _mm_add_pd(vec1_128, vec2_128);                                              \
        lasti += 1;                                                                            \
    }                                                                                          \
    found_nan = _mm256_movemask_pd(nan_flag) | _mm_movemask_pd(nan_flag128); \
    if(found_nan)                                                                              \
    {                                                                                          \
        max_value = NAN;                                                                       \
    }                                                                                          \
    else                                                                                       \
    {                                                                                          \
        /* Finding the maximum element in the vector */                                        \
        fla_get_max_xabs_element_real_avx2_max_reduce(d); /* Maximum value */                  \
        max_value = _mm_cvtsd_f64(vmax128);                                                    \
        max_value = sqrt(max_value);                                                           \
    }

#define fla_get_max_xabs_element_vector_complex_avx2_scaled(x, m, a, a_dim, scale_factor) \
    __m256d vec1, vec2, vec3, vec4;                                                       \
    __m256d vec5, vec6, vec7, vec8, scale_vec;                                            \
    __m128d vec1_128, vec2_128, scale_vec128;                                             \
    fla_get_max_xabs_element_vector_complex_avx2_define_extra_vars(x);                    \
    __m256d vmax = {0};                                                                   \
    __m128d vmax128 = {0};                                                                \
    doublereal max_value;                                                                 \
    integer i__, lasti;                                                                   \
    scale_vec = _mm256_set1_pd(scale_factor);                                             \
    scale_vec128 = _mm_set1_pd(scale_factor);                                             \
    max_value = 0.0f;                                                                     \
    lasti = m - (m % 8);                                                                  \
    /* Get maximum upto size lasti */                                                     \
    for(i__ = 1; i__ <= lasti; i__ += 8)                                                  \
    {                                                                                     \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec1, i__);               \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec2, i__ + 2);           \
        /* Scale values */                                                                \
        vec1 = _mm256_mul_pd(vec1, scale_vec);                                            \
        vec2 = _mm256_mul_pd(vec2, scale_vec);                                            \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                          \
        vec3 = _mm256_shuffle_pd(vec1, vec2, 0b0000);                                     \
        vec4 = _mm256_shuffle_pd(vec1, vec2, 0b1111);                                     \
        /* calculating r^2 and i^r */                                                     \
        vec3 = _mm256_mul_pd(vec3, vec3);                                                 \
        vec4 = _mm256_mul_pd(vec4, vec4);                                                 \
        /* vec3 = vec3 + vec4 */                                                          \
        vec3 = _mm256_add_pd(vec3, vec4);                                                 \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec5, i__ + 4);           \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec6, i__ + 6);           \
        /* Scale values */                                                                \
        vec5 = _mm256_mul_pd(vec5, scale_vec);                                            \
        vec6 = _mm256_mul_pd(vec6, scale_vec);                                            \
        /* vec7 = r = real(a[nx:2*nx]); vec8 = i = image(a[nx:2*nx])*/                    \
        vec7 = _mm256_shuffle_pd(vec5, vec6, 0b0000);                                     \
        vec8 = _mm256_shuffle_pd(vec5, vec6, 0b1111);                                     \
        /* calculating r^2 and i^r */                                                     \
        vec7 = _mm256_mul_pd(vec7, vec7);                                                 \
        vec8 = _mm256_mul_pd(vec8, vec8);                                                 \
        /* vec7 = vec7 + vec8 */                                                          \
        vec7 = _mm256_add_pd(vec7, vec8);                                                 \
        /* vec1 = max(vec3, vec7) */                                                      \
        vec1 = _mm256_max_pd(vec3, vec7);                                                 \
        vmax = _mm256_max_pd(vec1, vmax);                                                 \
    }                                                                                     \
    if((m - lasti) >= 4)                                                                  \
    {                                                                                     \
        /* Permute and find the maximum element*/                                         \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec1, lasti + 1);         \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec2, lasti + 3);         \
        /* Scale values */                                                                \
        vec1 = _mm256_mul_pd(vec1, scale_vec);                                            \
        vec2 = _mm256_mul_pd(vec2, scale_vec);                                            \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                          \
        vec3 = _mm256_shuffle_pd(vec1, vec2, 0b0000);                                     \
        vec4 = _mm256_shuffle_pd(vec1, vec2, 0b1111);                                     \
        /* calculating r^2 and i^r */                                                     \
        vec3 = _mm256_mul_pd(vec3, vec3);                                                 \
        vec4 = _mm256_mul_pd(vec4, vec4);                                                 \
        /* vec3 = vec3 + vec4 */                                                          \
        vec3 = _mm256_add_pd(vec3, vec4);                                                 \
        vmax = _mm256_max_pd(vec3, vmax);                                                 \
        lasti += 4;                                                                       \
    }                                                                                     \
    if((m - lasti) >= 2)                                                                  \
    {                                                                                     \
        fla_get_max_xabs_element_vector_complex_avx2_load256(x, vec1, lasti + 1);         \
        /* Scale values */                                                                \
        vec1 = _mm256_mul_pd(vec1, scale_vec);                                            \
        vec1 = _mm256_mul_pd(vec1, vec1);                                                 \
        vec2 = _mm256_permute_pd(vec1, 0b0101);                                           \
        vec1 = _mm256_add_pd(vec1, vec2);                                                 \
        vmax = _mm256_max_pd(vec1, vmax);                                                 \
        lasti += 2;                                                                       \
    }                                                                                     \
    if((m - lasti) == 1)                                                                  \
    {                                                                                     \
        fla_get_max_xabs_element_vector_complex_avx2_load128(x, vec1_128, m);             \
        /* Scale values */                                                                \
        vec1_128 = _mm_mul_pd(vec1_128, scale_vec128);                                    \
        vec1_128 = _mm_mul_pd(vec1_128, vec1_128);                                        \
        vec2_128 = _mm_permute_pd(vec1_128, 0b01);                                        \
        vmax128 = _mm_add_pd(vec1_128, vec2_128);                                         \
        lasti += 1;                                                                       \
    }                                                                                     \
    /* Finding the maximum element in the vector */                                       \
    fla_get_max_xabs_element_real_avx2_max_reduce(d); /* Maximum value */                 \
    max_value = _mm_cvtsd_f64(vmax128);                                                   \
    max_value = sqrt(max_value);

doublereal fla_get_max_dabs_element_vector_avx2(integer m, doublereal *a, integer a_dim)
{
    fla_get_max_xabs_element_vector_real_avx2(d, doublereal, m, a, a_dim, 4, 8, 2);
    return max_value;
}

real fla_get_max_sabs_element_vector_avx2(integer m, real *a, integer a_dim)
{
    fla_get_max_xabs_element_vector_real_avx2(s, real, m, a, a_dim, 8, 16, 4);
    return max_value;
}

doublereal fla_get_max_zabs_element_vector_avx2(integer m, doublecomplex *a, integer a_dim)
{
    /* Threshold amd scaling values taken from dlassq API
       which calculates sum of squares by performaing scaling
       to avoid overflow/underflow. */
    doublereal tsml = 1.4916681462400413E-154;
    doublereal tbig = 1.9979190722022350E+146;
    doublereal ssml = 4.4989137945431964E+161;
    doublereal sbig = 1.1113793747425387E-162;
    fla_get_max_xabs_element_vector_complex_avx2(d, m, a, a_dim);
    /* To avoid overflow if the max element is bigger/smaller than
       the threshold value then scale each element before
       calculating the absolute value.
       Also if the calculated value is NaN then no need to check */
    if(max_value > tbig && max_value == max_value)
    {
        fla_get_max_xabs_element_vector_complex_avx2_scaled(d, m, a, a_dim, sbig);
        return max_value / sbig;
    }
    else if(max_value < tsml && max_value == max_value)
    {
        fla_get_max_xabs_element_vector_complex_avx2_scaled(d, m, a, a_dim, ssml);
        return max_value / ssml;
    }
    return max_value;
}

real fla_get_max_cabs_element_vector_avx2(integer m, complex *a, integer a_dim)
{
    /* The real value is upscaled to doublereal during calculation of absolute
        value. Hence underflow/overflow is taken care of. */
    fla_get_max_xabs_element_vector_complex_avx2(s, m, a, a_dim);
    return max_value;
}
#endif
