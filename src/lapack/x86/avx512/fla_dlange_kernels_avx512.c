/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT
/* Find maxmimum absoute value of given vector using avx512 intrinsics*/
doublereal fla_get_max_abs_element_vector_avx512(integer m, doublereal *a, integer a_dim)
{
    __m512d vec1, vec2, vec3;
    __m512d vmax = {0};
    __m512d vd8_zero = _mm512_set1_pd(-0.0f);
    integer i__, lasti;
    doublereal max_value, temp;
    max_value = 0.0f;
    lasti = m - (m % 16);

    /* Find maxmimum absoute value loading 16 double elements
       in one iteration*/
    for(i__ = 1; i__ <= lasti; i__ += 16)
    {
        vec1 = _mm512_loadu_pd(&a[i__ + a_dim]);
        vec2 = _mm512_loadu_pd(&a[i__ + 8 + a_dim]);

        /* Find absolute values*/
        vec1 = _mm512_andnot_pd(vd8_zero, vec1);
        vec2 = _mm512_andnot_pd(vd8_zero, vec2);

        /* Find intermediate maximum values*/
        vec3 = _mm512_max_pd(vec1, vec2);
        vmax = _mm512_max_pd(vmax, vec3);
    }
    /* Maximum value at the end of iteration*/
    max_value = _mm512_reduce_max_pd(vmax);

    /* Find maximum from remainder loop*/
    for(i__ = (lasti + 1); i__ <= m; ++i__)
    {
        temp = f2c_abs(a[i__ + a_dim]);
        if(max_value < temp || temp != temp)
        {
            max_value = temp;
        }
    }
    return max_value;
}
#endif
