/******************************************************************************
 * * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

#define invoke_mm512_loadu_px(x, ...) _mm512_loadu_p##x(__VA_ARGS__)
#define invoke_mm512_max_px(x, ...) _mm512_max_p##x(__VA_ARGS__)
#define invoke_mm512_andnot_px(x, ...) _mm512_andnot_p##x(__VA_ARGS__)
#define invoke_mm512_set1_px(x, ...) _mm512_set1_p##x(__VA_ARGS__)
#define invoke_mm512_shuffle_px(x, ...) _mm512_shuffle_p##x(__VA_ARGS__)
#define invoke_mm512_mul_px(x, ...) _mm512_mul_p##x(__VA_ARGS__)
#define invoke_mm512_add_px(x, ...) _mm512_add_p##x(__VA_ARGS__)
#define invoke_mm512_reduce_max_px(x, ...) _mm512_reduce_max_p##x(__VA_ARGS__)
#define invoke_mm512_maskz_loadu_px(x, ...) _mm512_maskz_loadu_p##x(__VA_ARGS__)
#define invoke_mm512_permute_px(x, ...) _mm512_permute_p##x(__VA_ARGS__)

#define invoke_sqrts(data) sqrtf(data)
#define invoke_sqrtd(data) sqrt(data)
#define invoke_sqrtx(x, data) invoke_sqrt##x(data)

#define type_m512s __m512
#define type_m512d __m512d
#define type_mm512x(x) type_m512##x

#define fla_get_max_xabs_element_real_avx512(x, realtype, m, a, a_dim, nx, hnx, stepsize) \
    type_mm512x(x) vec1, vec2, vec3;                                                      \
    type_mm512x(x) vmax = {0};                                                            \
    type_mm512x(x) vd8_nzero = invoke_mm512_set1_px(x, -0.0f);                            \
    integer i__, lasti, remain;                                                           \
    realtype max_value;                                                                   \
    max_value = 0.0f;                                                                     \
    lasti = m - (m % stepsize);                                                           \
    /* Get maximum upto size lasti */                                                     \
    for(i__ = 1; i__ <= lasti; i__ += stepsize)                                           \
    {                                                                                     \
        vec1 = invoke_mm512_loadu_px(x, &a[i__ + a_dim]);                                 \
        vec2 = invoke_mm512_loadu_px(x, &a[i__ + nx + a_dim]);                            \
        /* Find absolute values */                                                        \
        vec1 = invoke_mm512_andnot_px(x, vd8_nzero, vec1);                                \
        vec2 = invoke_mm512_andnot_px(x, vd8_nzero, vec2);                                \
        /* Find intermediate maximum values */                                            \
        vec3 = invoke_mm512_max_px(x, vec1, vec2);                                        \
        vmax = invoke_mm512_max_px(x, vec3, vmax);                                        \
    }                                                                                     \
    /* reduce more nx elements*/                                                          \
    if((m - lasti) >= nx)                                                                 \
    {                                                                                     \
        vec1 = invoke_mm512_loadu_px(x, &a[lasti + 1 + a_dim]);                           \
        vec1 = invoke_mm512_andnot_px(x, vd8_nzero, vec1);                                \
        vmax = invoke_mm512_max_px(x, vec1, vmax);                                        \
        lasti += nx;                                                                      \
    }                                                                                     \
    if(lasti < m)                                                                         \
    {                                                                                     \
        remain = m - lasti;                                                               \
        vec1 = invoke_mm512_maskz_loadu_px(x, (1 << remain) - 1, &a[lasti + 1 + a_dim]);  \
        vec1 = invoke_mm512_andnot_px(x, vd8_nzero, vec1);                                \
        vmax = invoke_mm512_max_px(x, vec1, vmax);                                        \
    }                                                                                     \
    max_value = invoke_mm512_reduce_max_px(x, vmax);

#define fla_get_max_sabs_element_complex_mm512_shuffle_flag_real 0b10001000
#define fla_get_max_sabs_element_complex_mm512_shuffle_flag_imag 0b11011101

#define fla_get_max_dabs_element_complex_mm512_shuffle_flag_real 0b00000000
#define fla_get_max_dabs_element_complex_mm512_shuffle_flag_imag 0b11111111

#define fla_get_max_sabs_element_complex_mm512_shuffle_ri_flag 0b10110001
#define fla_get_max_dabs_element_complex_mm512_shuffle_ri_flag 0b01010101

#define fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(x) \
    fla_get_max_##x##abs_element_complex_mm512_shuffle_flag_real

#define fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(x) \
    fla_get_max_##x##abs_element_complex_mm512_shuffle_flag_imag

#define fla_get_max_xabs_element_complex_mm512_shuffle_ri_flag(x) \
    fla_get_max_##x##abs_element_complex_mm512_shuffle_ri_flag

#define fla_get_max_sabs_element_vector_complex_avx512_extra_vars \
    __m256 vec1_256s, vec2_256s, vec5_256s, vec6_256s;

#define fla_get_max_dabs_element_vector_complex_avx512_extra_vars

#define fla_get_max_xabs_element_vector_complex_avx512_define_extra_vars(x) \
    fla_get_max_##x##abs_element_vector_complex_avx512_extra_vars

#define fla_get_max_sabs_element_vector_complex_avx2_load512(target_vec, idx) \
    target_vec##_256s = _mm256_loadu_ps((real *)&a[(idx) + a_dim]);                  \
    target_vec = _mm512_cvtps_pd(target_vec##_256s);

#define fla_get_max_dabs_element_vector_complex_avx2_load512(target_vec, idx) \
    target_vec = _mm512_loadu_pd((doublereal *)&a[(idx) + a_dim]);

#define fla_get_max_xabs_element_vector_complex_avx2_load512(x, target_vec, idx) \
    fla_get_max_##x##abs_element_vector_complex_avx2_load512(target_vec, idx)

#define fla_get_max_sabs_element_vector_complex_avx2_load512_masked(target_vec, idx, mask) \
    target_vec##_256s                                                                             \
        = _mm512_castps512_ps256(_mm512_maskz_loadu_ps(mask, (real *)&a[(idx) + a_dim]));         \
    target_vec = _mm512_cvtps_pd(target_vec##_256s);

#define fla_get_max_dabs_element_vector_complex_avx2_load512_masked(target_vec, idx, mask) \
    target_vec = _mm512_maskz_loadu_pd(mask, (doublereal *)&a[(idx) + a_dim]);

#define fla_get_max_xabs_element_vector_complex_avx2_load512_masked(x, target_vec, idx, \
                                                                           mask)               \
    fla_get_max_##x##abs_element_vector_complex_avx2_load512_masked(target_vec, idx, mask)


#define fla_get_max_xabs_element_vector_complex_avx512(x, m, a, a_dim, nx, hnx,      \
                                                       stepsize)                               \
    type_mm512x(d) vec1, vec2, vec3, vec4;                                                     \
    type_mm512x(d) vec5, vec6, vec7, vec8;                                                     \
    fla_get_max_xabs_element_vector_complex_avx512_define_extra_vars(x);                       \
    type_mm512x(d) vmax = {0};                                                                 \
    doublereal max_value;                                                                        \
    max_value = 0.0;                                                                          \
    integer i__, lasti, remainv1, remainv2;                                                    \
    lasti = m - (m % 16);                                                                \
    /* Get maximum upto size lasti */                                                          \
    for(i__ = 1; i__ <= lasti; i__ += stepsize)                                                  \
    {                                                                                            \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec1, i__);               \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec2, i__ + hnx);         \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                                 \
        vec3 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));         \
        vec4 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));         \
        /* calculating r^2 and i^r */                                                            \
        vec3 = invoke_mm512_mul_px(d, vec3, vec3);                                               \
        vec4 = invoke_mm512_mul_px(d, vec4, vec4);                                               \
        /* vec3 = vec3 + vec4 */                                                                 \
        vec3 = invoke_mm512_add_px(d, vec3, vec4);                                               \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec5, i__ + nx);          \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec6, i__ + nx + hnx);    \
        /* vec7 = r = real(a[nx:2*nx]); vec8 = i = image(a[nx:2*nx])*/                           \
        vec7 = invoke_mm512_shuffle_px(                                                          \
            d, vec5, vec6, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));         \
        vec8 = invoke_mm512_shuffle_px(                                                          \
            d, vec5, vec6, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));         \
        /* calculating r^2 and i^r */                                                            \
        vec7 = invoke_mm512_mul_px(d, vec7, vec7);                                               \
        vec8 = invoke_mm512_mul_px(d, vec8, vec8);                                               \
        /* vec7 = vec7 + vec8 */                                                                 \
        vec7 = invoke_mm512_add_px(d, vec7, vec8);                                               \
        /* vec1 = max(vec3, vec7) */                                                             \
        vec1 = invoke_mm512_max_px(d, vec3, vec7);                                               \
        vmax = invoke_mm512_max_px(d, vec1, vmax);                                               \
    }                                                                                            \
    if((m - lasti) >= nx)                                                                        \
    {                                                                                            \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec1, lasti + 1);         \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec2, lasti + 1 + hnx);   \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                                 \
        vec3 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));         \
        vec4 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));         \
        /* calculating r^2 and i^r */                                                            \
        vec3 = invoke_mm512_mul_px(d, vec3, vec3);                                               \
        vec4 = invoke_mm512_mul_px(d, vec4, vec4);                                               \
        /* vec3 = vec3 + vec4 */                                                                 \
        vec3 = invoke_mm512_add_px(d, vec3, vec4);                                               \
        vmax = invoke_mm512_max_px(d, vec3, vmax);                                               \
        lasti += nx;                                                                             \
    }                                                                                            \
    if(lasti < m)                                                                                \
    {                                                                                            \
        remainv1 = fla_min(m - lasti, hnx);                                                      \
        remainv2 = m - lasti - remainv1;                                                         \
        remainv1 <<= 1;                                                                          \
        remainv2 <<= 1;                                                                          \
        if(remainv2 != 0)                                                                        \
        {                                                                                        \
            fla_get_max_xabs_element_vector_complex_avx2_load512_masked(                  \
                x, vec1, lasti + 1, ((1 << remainv1) - 1));                                      \
            fla_get_max_xabs_element_vector_complex_avx2_load512_masked(                  \
                x, vec2, lasti + 1 + hnx, ((1 << remainv2) - 1));                                \
            /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                             \
            vec3 = invoke_mm512_shuffle_px(                                                      \
                d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));     \
            vec4 = invoke_mm512_shuffle_px(                                                      \
                d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));     \
            /* calculating r^2 and i^r */                                                        \
            vec3 = invoke_mm512_mul_px(d, vec3, vec3);                                           \
            vec4 = invoke_mm512_mul_px(d, vec4, vec4);                                           \
            /* vec3 = vec3 + vec4 */                                                             \
            vec3 = invoke_mm512_add_px(d, vec3, vec4);                                           \
            vmax = invoke_mm512_max_px(d, vec3, vmax);                                           \
        }                                                                                        \
        else                                                                                     \
        {                                                                                        \
            fla_get_max_xabs_element_vector_complex_avx2_load512_masked(                  \
                x, vec1, lasti + 1, ((1 << remainv1) - 1));                                      \
            vec1 = invoke_mm512_mul_px(d, vec1, vec1);                                           \
            vec2 = invoke_mm512_permute_px(                                                      \
                d, vec1, fla_get_max_xabs_element_complex_mm512_shuffle_ri_flag(d));             \
            vec1 = invoke_mm512_add_px(d, vec1, vec2);                                           \
            vmax = invoke_mm512_max_px(d, vec1, vmax);                                           \
        }                                                                                        \
    }                                                                                            \
    max_value = invoke_mm512_reduce_max_px(d, vmax);                                             \
    max_value = invoke_sqrtx(d, max_value);


#define fla_get_max_xabs_element_vector_complex_avx512_scaled(x, m, a, a_dim, nx, hnx, \
                                                              stepsize, scale_factor)            \
    type_mm512x(d) vec1, vec2, vec3, vec4;                                                       \
    type_mm512x(d) vec5, vec6, vec7, vec8;                                                       \
    type_mm512x(d) vmax = {0}, scale_vec;                                                        \
    fla_get_max_xabs_element_vector_complex_avx512_define_extra_vars(x);                  \
    doublereal max_value;                                                                        \
    integer i__, lasti, remainv1, remainv2;                                                      \
    max_value = 0.0f;                                                                            \
    lasti = m - (m % stepsize);                                                                  \
    scale_vec = invoke_mm512_set1_px(d, scale_factor);                                           \
    /* Get maximum upto size lasti */                                                            \
    for(i__ = 1; i__ <= lasti; i__ += stepsize)                                                  \
    {                                                                                            \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec1, i__);               \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec2, i__ + hnx);         \
        /* Scale values */                                                                       \
        vec1 = invoke_mm512_mul_px(d, vec1, scale_vec);                                          \
        vec2 = invoke_mm512_mul_px(d, vec2, scale_vec);                                          \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                                 \
        vec3 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));         \
        vec4 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));         \
        /* calculating r^2 and i^r */                                                            \
        vec3 = invoke_mm512_mul_px(d, vec3, vec3);                                               \
        vec4 = invoke_mm512_mul_px(d, vec4, vec4);                                               \
        /* vec3 = vec3 + vec4 */                                                                 \
        vec3 = invoke_mm512_add_px(d, vec3, vec4);                                               \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec5, i__ + nx);          \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec6, i__ + nx + hnx);    \
        /* Scale values */                                                                       \
        vec5 = invoke_mm512_mul_px(d, vec5, scale_vec);                                          \
        vec6 = invoke_mm512_mul_px(d, vec6, scale_vec);                                          \
        /* vec7 = r = real(a[nx:2*nx]); vec8 = i = image(a[nx:2*nx])*/                           \
        vec7 = invoke_mm512_shuffle_px(                                                          \
            d, vec5, vec6, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));         \
        vec8 = invoke_mm512_shuffle_px(                                                          \
            d, vec5, vec6, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));         \
        /* calculating r^2 and i^r */                                                            \
        vec7 = invoke_mm512_mul_px(d, vec7, vec7);                                               \
        vec8 = invoke_mm512_mul_px(d, vec8, vec8);                                               \
        /* vec7 = vec7 + vec8 */                                                                 \
        vec7 = invoke_mm512_add_px(d, vec7, vec8);                                               \
        /* vec1 = max(vec3, vec7) */                                                             \
        vec1 = invoke_mm512_max_px(d, vec3, vec7);                                               \
        vmax = invoke_mm512_max_px(d, vec1, vmax);                                               \
    }                                                                                            \
    if((m - lasti) >= nx)                                                                        \
    {                                                                                            \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec1, lasti + 1);         \
        fla_get_max_xabs_element_vector_complex_avx2_load512(x, vec2, lasti + 1 + hnx);   \
        /* Scale values */                                                                       \
        vec1 = invoke_mm512_mul_px(d, vec1, scale_vec);                                          \
        vec2 = invoke_mm512_mul_px(d, vec2, scale_vec);                                          \
        /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                                 \
        vec3 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));         \
        vec4 = invoke_mm512_shuffle_px(                                                          \
            d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));         \
        /* calculating r^2 and i^r */                                                            \
        vec3 = invoke_mm512_mul_px(d, vec3, vec3);                                               \
        vec4 = invoke_mm512_mul_px(d, vec4, vec4);                                               \
        /* vec3 = vec3 + vec4 */                                                                 \
        vec3 = invoke_mm512_add_px(d, vec3, vec4);                                               \
        vmax = invoke_mm512_max_px(d, vec3, vmax);                                               \
        lasti += nx;                                                                             \
    }                                                                                            \
    if(lasti < m)                                                                                \
    {                                                                                            \
        remainv1 = fla_min(m - lasti, hnx);                                                      \
        remainv2 = m - lasti - remainv1;                                                         \
        remainv1 <<= 1;                                                                          \
        remainv2 <<= 1;                                                                          \
        if(remainv2 != 0)                                                                        \
        {                                                                                        \
            fla_get_max_xabs_element_vector_complex_avx2_load512_masked(                  \
                x, vec1, lasti + 1, ((1 << remainv1) - 1));                                      \
            fla_get_max_xabs_element_vector_complex_avx2_load512_masked(                  \
                x, vec2, lasti + 1 + hnx, ((1 << remainv2) - 1));                                \
            /* Scale values */                                                                   \
            vec1 = invoke_mm512_mul_px(d, vec1, scale_vec);                                      \
            vec2 = invoke_mm512_mul_px(d, vec2, scale_vec);                                      \
            /* vec3 = r = real(a[0:nx]); vec4 = i = image(a[0:nx])*/                             \
            vec3 = invoke_mm512_shuffle_px(                                                      \
                d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_real(d));     \
            vec4 = invoke_mm512_shuffle_px(                                                      \
                d, vec1, vec2, fla_get_max_xabs_element_complex_mm512_shuffle_flag_imag(d));     \
            /* calculating r^2 and i^r */                                                        \
            vec3 = invoke_mm512_mul_px(d, vec3, vec3);                                           \
            vec4 = invoke_mm512_mul_px(d, vec4, vec4);                                           \
            /* vec3 = vec3 + vec4 */                                                             \
            vec3 = invoke_mm512_add_px(d, vec3, vec4);                                           \
            vmax = invoke_mm512_max_px(d, vec3, vmax);                                           \
        }                                                                                        \
        else                                                                                     \
        {                                                                                        \
            fla_get_max_xabs_element_vector_complex_avx2_load512_masked(                  \
                x, vec1, lasti + 1, ((1 << remainv1) - 1));                                      \
            /* Scale values */                                                                   \
            vec1 = invoke_mm512_mul_px(d, vec1, scale_vec);                                      \
            vec1 = invoke_mm512_mul_px(d, vec1, vec1);                                           \
            vec2 = invoke_mm512_permute_px(                                                      \
                d, vec1, fla_get_max_xabs_element_complex_mm512_shuffle_ri_flag(d));             \
            vec1 = invoke_mm512_add_px(d, vec1, vec2);                                           \
            vmax = invoke_mm512_max_px(d, vec1, vmax);                                           \
        }                                                                                        \
    }                                                                                            \
    max_value = invoke_mm512_reduce_max_px(d, vmax);                                             \
    max_value = invoke_sqrtx(d, max_value);

/* Find maxmimum absoute value of given doublereal vector using avx512 intrinsics*/
doublereal fla_get_max_dabs_element_vector_avx512(integer m, doublereal *a, integer a_dim)
{
    fla_get_max_xabs_element_real_avx512(d, doublereal, m, a, a_dim, 8, 4, 16);
    return max_value;
}

/* Find maxmimum absoute value of given real vector using avx512 intrinsics*/
real fla_get_max_sabs_element_vector_avx512(integer m, real *a, integer a_dim)
{
    fla_get_max_xabs_element_real_avx512(s, real, m, a, a_dim, 16, 8, 32);
    return max_value;
}

/* Find maxmimum absoute value of given doublecomplex vector using avx512 intrinsics*/
doublereal fla_get_max_zabs_element_vector_avx512(integer m, doublecomplex *a, integer a_dim)
{
    /* Threshold amd scaling values taken from dlassq API
       which calculates sum of squares by performaing scaling
       to avoid overflow/underflow. */
    doublereal tsml = 1.4916681462400413E-154;
    doublereal tbig = 1.9979190722022350E+146;
    doublereal ssml = 4.4989137945431964E+161;
    doublereal sbig = 1.1113793747425387E-162;
    fla_get_max_xabs_element_vector_complex_avx512(d, m, a, a_dim, 8, 4, 16);
    /* To avoid overflow if the max element is bigger/smaller than
       the threshold value then scale each element before
       calculating the absolute value. */
    if(max_value > tbig)
    {
        fla_get_max_xabs_element_vector_complex_avx512_scaled(d, m, a, a_dim, 8, 4, 16,
                                                              sbig);
        return max_value / sbig;
    }
    else if(max_value < tsml)
    {
        fla_get_max_xabs_element_vector_complex_avx512_scaled(d, m, a, a_dim, 8, 4, 16,
                                                              ssml);
        return max_value / ssml;
    }
    return max_value;
}

/* Find maxmimum absoute value of given complex vector using avx512 intrinsics*/
real fla_get_max_cabs_element_vector_avx512(integer m, complex *a, integer a_dim)
{
    /* The real value is upscaled to doublereal during calculation of absolute
       value. Hence underflow/overflow is taken care of. */
    fla_get_max_xabs_element_vector_complex_avx512(s, m, a, a_dim, 8, 4, 16);
    return max_value;
}

#endif
