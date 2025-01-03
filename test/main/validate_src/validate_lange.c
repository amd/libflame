/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_lange.c
 *  @brief Defines validate function of LANGE() to use in test suite.
 *  */

#include "test_common.h"

/* Used for scaling during sum of squares calculation */
#define float_sbig 1.32348898E-23f
#define float_ssml 3.77789319E+22f
#define float_tsml 1.08420217E-19f
#define float_tbig 4.50359963E+15f
#define double_tsml 1.4916681462400413E-154
#define double_tbig 1.9979190722022350E+146
#define double_ssml 4.4989137945431964E+161
#define double_sbig 1.1113793747425387E-162

#define x_sbig(x) x##_sbig
#define x_ssml(x) x##_ssml
#define x_tsml(x) x##_tsml
#define x_tbig(x) x##_tbig

#ifdef _WIN32
#define GET_ABS_VAL(datatype, A, i, j, lda, result)                                               \
    {                                                                                             \
        switch(datatype)                                                                          \
        {                                                                                         \
            case FLOAT:                                                                           \
            {                                                                                     \
                result = FLA_FABS(((float *)A)[i + j * lda]);                                     \
                break;                                                                            \
            }                                                                                     \
            case DOUBLE:                                                                          \
            {                                                                                     \
                result = FLA_FABS(((double *)A)[i + j * lda]);                                    \
                break;                                                                            \
            }                                                                                     \
            case COMPLEX:                                                                         \
            {                                                                                     \
                _Fcomplex z_ = {((complex *)A)[i + j * lda].r, ((complex *)A)[i + j * lda].i};    \
                result = cabsf(z_);                                                               \
                break;                                                                            \
            }                                                                                     \
            case DOUBLE_COMPLEX:                                                                  \
            {                                                                                     \
                _Dcomplex z_                                                                      \
                    = {((doublecomplex *)A)[i + j * lda].r, ((doublecomplex *)A)[i + j * lda].i}; \
                result = cabs(z_);                                                                \
                break;                                                                            \
            }                                                                                     \
        }                                                                                         \
    }
#else
#define GET_ABS_VAL(datatype, A, i, j, lda, result)                                               \
    {                                                                                             \
        switch(datatype)                                                                          \
        {                                                                                         \
            case FLOAT:                                                                           \
                result = FLA_FABS(((float *)A)[i + j * lda]);                                     \
                break;                                                                            \
            case DOUBLE:                                                                          \
                result = FLA_FABS(((double *)A)[i + j * lda]);                                    \
                break;                                                                            \
            case COMPLEX:                                                                         \
                result = cabs(((complex *)A)[i + j * lda].r + I * ((complex *)A)[i + j * lda].i); \
                break;                                                                            \
            case DOUBLE_COMPLEX:                                                                  \
                result = cabs(((doublecomplex *)A)[i + j * lda].r                                 \
                              + I * ((doublecomplex *)A)[i + j * lda].i);                         \
                break;                                                                            \
        }                                                                                         \
    }
#endif

#define GET_REAL_ABS_VAL(datatype, A, i, j, lda, result)                \
    {                                                                   \
        switch(datatype)                                                \
        {                                                               \
            case FLOAT:                                                 \
                result = FLA_FABS(((float *)A)[i + j * lda]);           \
                break;                                                  \
            case DOUBLE:                                                \
                result = FLA_FABS(((double *)A)[i + j * lda]);          \
                break;                                                  \
            case COMPLEX:                                               \
                result = FLA_FABS(((complex *)A)[i + j * lda].r);       \
                break;                                                  \
            case DOUBLE_COMPLEX:                                        \
                result = FLA_FABS(((doublecomplex *)A)[i + j * lda].r); \
                break;                                                  \
        }                                                               \
    }

#define GET_IMAG_ABS_VAL(datatype, A, i, j, lda, result)                \
    {                                                                   \
        switch(datatype)                                                \
        {                                                               \
            case COMPLEX:                                               \
                result = FLA_FABS(((complex *)A)[i + j * lda].i);       \
                break;                                                  \
            case DOUBLE_COMPLEX:                                        \
                result = FLA_FABS(((doublecomplex *)A)[i + j * lda].i); \
                break;                                                  \
            default:                                                    \
                result = 0.0;                                           \
                break;                                                  \
        }                                                               \
    }

#define GET_SQR_VAL(datatype, A, i, j, lda, result)                                   \
    {                                                                                 \
        switch(datatype)                                                              \
        {                                                                             \
            case FLOAT:                                                               \
            {                                                                         \
                result = (((float *)A)[i + j * lda]) * (((float *)A)[i + j * lda]);   \
                break;                                                                \
            }                                                                         \
            case DOUBLE:                                                              \
            {                                                                         \
                result = (((double *)A)[i + j * lda]) * (((double *)A)[i + j * lda]); \
                break;                                                                \
            }                                                                         \
            case COMPLEX:                                                             \
            {                                                                         \
                complex *C = (complex *)A;                                            \
                result = (C[i + j * lda].r * C[i + j * lda].r)                        \
                         + (C[i + j * lda].i * C[i + j * lda].i);                     \
                break;                                                                \
            }                                                                         \
            case DOUBLE_COMPLEX:                                                      \
            {                                                                         \
                doublecomplex *C = (doublecomplex *)A;                                \
                result = (C[i + j * lda].r * C[i + j * lda].r)                        \
                         + (C[i + j * lda].i * C[i + j * lda].i);                     \
                break;                                                                \
            }                                                                         \
        }                                                                             \
    }

#define GET_MNORM(realtype, datatype, A, m, n, lda, resultp)  \
    {                                                         \
        integer i, j;                                         \
        realtype max_val = 0.0;                               \
        realtype abs_val;                                     \
        for(j = 0; j < n; j++)                                \
        {                                                     \
            for(i = 0; i < m; i++)                            \
            {                                                 \
                GET_ABS_VAL(datatype, A, i, j, lda, abs_val); \
                max_val = fla_max(max_val, abs_val);          \
            }                                                 \
        }                                                     \
        *(realtype *)resultp = max_val;                       \
    }

#define GET_1NORM(realtype, datatype, A, m, n, lda, resultp)  \
    {                                                         \
        integer i, j;                                         \
        realtype max_val = 0.0, col_sum = 0.0;                \
        double abs_val;                                       \
        for(j = 0; j < n; j++)                                \
        {                                                     \
            col_sum = 0.0;                                    \
            for(i = 0; i < m; i++)                            \
            {                                                 \
                GET_ABS_VAL(datatype, A, i, j, lda, abs_val); \
                col_sum += abs_val;                           \
            }                                                 \
            max_val = fla_max(max_val, col_sum);              \
        }                                                     \
        *(realtype *)resultp = max_val;                       \
    }

#define GET_INORM(realtype, datatype, A, m, n, lda, resultp)          \
    {                                                                 \
        integer i, j;                                                 \
        realtype max_val = 0.0;                                       \
        double abs_val;                                               \
        realtype *row_sums;                                           \
        create_vector(get_realtype(datatype), (void **)&row_sums, m); \
        reset_vector(get_realtype(datatype), row_sums, m, 1);         \
        for(j = 0; j < n; j++)                                        \
        {                                                             \
            for(i = 0; i < m; i++)                                    \
            {                                                         \
                GET_ABS_VAL(datatype, A, i, j, lda, abs_val);         \
                row_sums[i] += abs_val;                               \
            }                                                         \
        }                                                             \
        for(i = 0; i < m; i++)                                        \
        {                                                             \
            max_val = fla_max(max_val, row_sums[i]);                  \
        }                                                             \
        free_vector(row_sums);                                        \
        *(realtype *)resultp = max_val;                               \
    }

#define GET_FNORM(realtype, datatype, A, m, n, lda, resultp)            \
    {                                                                   \
        integer i, j;                                                   \
        integer notbig = 1;                                             \
        realtype sqrsum = 0., scl = 1.;                                 \
        realtype abig, amed, asml, abs_val, t__;                        \
        for(j = 0; j < n; j++)                                          \
        {                                                               \
            abig = 0.;                                                  \
            amed = 0.;                                                  \
            asml = 0.;                                                  \
            for(i = 0; i < m; i++)                                      \
            {                                                           \
                GET_REAL_ABS_VAL(datatype, A, i, j, lda, abs_val);      \
                if(abs_val > x_tbig(realtype))                          \
                {                                                       \
                    t__ = abs_val * x_sbig(realtype);                   \
                    abig += t__ * t__;                                  \
                    notbig = 0;                                         \
                }                                                       \
                else if(abs_val < x_tsml(realtype))                     \
                {                                                       \
                    if(notbig)                                          \
                    {                                                   \
                        t__ = abs_val * x_ssml(realtype);               \
                        asml += t__ * t__;                              \
                    }                                                   \
                }                                                       \
                else                                                    \
                {                                                       \
                    amed += abs_val * abs_val;                          \
                }                                                       \
                if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)   \
                {                                                       \
                    GET_IMAG_ABS_VAL(datatype, A, i, j, lda, abs_val);  \
                    if(abs_val > x_tbig(realtype))                      \
                    {                                                   \
                        t__ = abs_val * x_sbig(realtype);               \
                        abig += t__ * t__;                              \
                        notbig = 0;                                     \
                    }                                                   \
                    else if(abs_val < x_tsml(realtype))                 \
                    {                                                   \
                        if(notbig)                                      \
                        {                                               \
                            t__ = abs_val * x_ssml(realtype);           \
                            asml += t__ * t__;                          \
                        }                                               \
                    }                                                   \
                    else                                                \
                    {                                                   \
                        amed += abs_val * abs_val;                      \
                    }                                                   \
                }                                                       \
            }                                                           \
            if(sqrsum > 0)                                              \
            {                                                           \
                abs_val = scl * sqrt(sqrsum);                           \
                if(abs_val > x_tbig(realtype))                          \
                {                                                       \
                    t__ = scl * x_sbig(realtype);                       \
                    abig += (t__ * t__) * sqrsum;                       \
                    notbig = 0;                                         \
                }                                                       \
                else if(abs_val < x_tsml(realtype))                     \
                {                                                       \
                    if(notbig)                                          \
                    {                                                   \
                        t__ = scl * x_ssml(realtype);                   \
                        asml += (t__ * t__) * sqrsum;                   \
                    }                                                   \
                }                                                       \
                else                                                    \
                {                                                       \
                    amed += (scl * scl) * sqrsum;                       \
                }                                                       \
            }                                                           \
            if(abig > 0.)                                               \
            {                                                           \
                if(amed > 0. || amed != amed)                           \
                {                                                       \
                    abig += amed * x_sbig(realtype) * x_sbig(realtype); \
                }                                                       \
                scl = (realtype)1. / x_sbig(realtype);                  \
                sqrsum = abig;                                          \
            }                                                           \
            else if(asml > 0.)                                          \
            {                                                           \
                if(amed > 0. || amed != amed)                           \
                {                                                       \
                    amed = sqrt(amed);                                  \
                    asml = sqrt(asml) / x_ssml(realtype);               \
                    realtype ymax = fla_max(amed, asml);                \
                    realtype ymin = fla_min(amed, asml);                \
                    realtype t1 = ymin / ymax;                          \
                    scl = 1.;                                           \
                    sqrsum = ymax * ymax * (t1 * t1 + 1.);              \
                }                                                       \
                else                                                    \
                {                                                       \
                    sqrsum = asml;                                      \
                    scl = (realtype)1. / x_ssml(realtype);              \
                }                                                       \
            }                                                           \
            else                                                        \
            {                                                           \
                sqrsum = amed;                                          \
                scl = 1.;                                               \
            }                                                           \
        }                                                               \
        *(realtype *)resultp = scl * sqrt(sqrsum);                      \
    }

#define CALL_NORM(datatype, NORM_MACRO, ...)                     \
    {                                                            \
        switch(datatype)                                         \
        {                                                        \
            case FLOAT:                                          \
            {                                                    \
                NORM_MACRO(float, FLOAT, __VA_ARGS__);           \
                break;                                           \
            }                                                    \
            case DOUBLE:                                         \
            {                                                    \
                NORM_MACRO(double, DOUBLE, __VA_ARGS__);         \
                break;                                           \
            }                                                    \
            case COMPLEX:                                        \
            {                                                    \
                NORM_MACRO(float, COMPLEX, __VA_ARGS__);         \
                break;                                           \
            }                                                    \
            case DOUBLE_COMPLEX:                                 \
            {                                                    \
                NORM_MACRO(double, DOUBLE_COMPLEX, __VA_ARGS__); \
                break;                                           \
            }                                                    \
        }                                                        \
    }

void validate_lange(integer datatype, char norm_type, integer m, integer n, integer lda, void *A,
                    void *result, double *residual)
{
    void *calculated_value;

    create_vector(get_realtype(datatype), &calculated_value, 1);

    /* Calculate the target norm value */
    switch(norm_type)
    {
        case 'M':
        {
            CALL_NORM(datatype, GET_MNORM, A, m, n, lda, calculated_value);
            break;
        }
        case '1':
        {
            CALL_NORM(datatype, GET_1NORM, A, m, n, lda, calculated_value);
            break;
        }
        case 'I':
        {
            CALL_NORM(datatype, GET_INORM, A, m, n, lda, calculated_value);
            break;
        }
        case 'F':
        {
            CALL_NORM(datatype, GET_FNORM, A, m, n, lda, calculated_value);
            break;
        }
    }

    /* Calculate the residual */
    switch(get_realtype(datatype))
    {
        case FLOAT:
        {
            double res_value = *(float *)result - *(float *)calculated_value;
            *residual = FLA_FABS(res_value);
            float eps = slamch_("P");
            if(norm_type == 'F' && *residual > (m * n * eps))
            {
                *residual = DBL_MAX;
            }
            else if(norm_type != 'F' && *residual > eps)
            {
                *residual = DBL_MAX;
            }
            break;
        }
        case DOUBLE:
        {
            double res_value = *(double *)result - *(double *)calculated_value;
            *residual = FLA_FABS(res_value);
            double eps = dlamch_("P");
            if(norm_type == 'F' && *residual > (m * n * eps))
            {
                *residual = DBL_MAX;
            }
            else if(norm_type != 'F' && *residual > eps)
            {
                *residual = DBL_MAX;
            }
            break;
        }
    }

    free_vector(calculated_value);
}