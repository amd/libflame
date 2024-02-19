/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
void fla_dlarf_left_apply_incv1_avx2(integer m, integer n, doublereal *r, integer ldr,
                                     doublereal *v, doublereal ntau, doublereal *work)
{
    integer acols, arows;
    integer k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp, vd2_vj1;
    __m256d vd4_dtmp, vd4_inp, vd4_vj;
    __m128d vd2_ltmp, vd2_htmp;

    /* Part 2: Apply the Householder rotation              */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */

    arows = m;
    acols = n;
    vd2_ntau = _mm_set1_pd(ntau);

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp = _mm_setzero_pd();
        vd4_dtmp = _mm256_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&r[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp = _mm256_fmadd_pd(vd4_inp, vd4_vj, vd4_dtmp);
        }
        if(k < arows)
        {
            /* load column elements of A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&r[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp = _mm_fmadd_pd(vd2_inp, vd2_vj1, vd2_dtmp);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&r[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp = _mm_fmadd_pd(vd2_inp, vd2_vj1, vd2_dtmp);
        }
        /* Horizontal add of dtmp */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);
        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        /* alternate for above 2 instructions which do not  */
        /* compile for older gcc versions (7 and below).    */
        /* Both will be same in terms of latency though     */
        /* vd4_dtmp = _mm256_set_m128d(vd2_dtmp, vd2_dtmp); */

        /* Compute c_A + tmp * v */
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&r[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_inp = _mm256_fmadd_pd(vd4_dtmp, vd4_vj, vd4_inp);
            _mm256_storeu_pd((doublereal *)&r[k + j * ldr], vd4_inp);
        }
        if(k < arows)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&r[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_inp = _mm_fmadd_pd(vd2_dtmp, vd2_vj1, vd2_inp);
            _mm_storeu_pd((doublereal *)&r[k + j * ldr], vd2_inp);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&r[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_inp = _mm_fmadd_pd(vd2_dtmp, vd2_vj1, vd2_inp);
            _mm_storel_pd((doublereal *)&r[k + j * ldr], vd2_inp);
        }
    }
}
#endif