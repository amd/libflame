/******************************************************************************
 * * Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
__attribute__((aligned(512))) void fla_dlarf_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n,
                                                                   doublereal *a_buff,
                                                                   aocl_int64_t ldr, doublereal *v,
                                                                   doublereal ntau,
                                                                   doublereal *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp, vd2_vj1, vd2_dtmp2;
    __m256d vd4_dtmp, vd4_inp, vd4_vj, vd4_dtmp2;
    __m128d vd2_ltmp, vd2_htmp;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    arows = m;
    acols = n;
    vd2_ntau = _mm_set1_pd(ntau);
    --v;
    a_buff -= ldr + 1;
    --work;

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp = _mm_setzero_pd();
        vd4_dtmp = _mm256_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);
        }
        if(k < arows)
        {
            /* load column elements of A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
        }
        /* Horizontal add of dtmp */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);
        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        /* Store the result in work */
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
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
        }
        if(k < arows)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storel_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }
    }
}

/* v(1) is not stored explicitly and is assumed to be 1 */
__attribute__((aligned(512))) void fla_dlarf1f_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n,
                                                                     doublereal *a_buff,
                                                                     aocl_int64_t ldr,
                                                                     doublereal *v, doublereal ntau,
                                                                     doublereal *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp, vd2_vj1, vd2_dtmp2, vd2_dtmp3;
    __m256d vd4_dtmp, vd4_inp, vd4_vj, vd4_dtmp2;
    __m128d vd2_ltmp, vd2_htmp;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    arows = m;
    acols = n;
    vd2_ntau = _mm_set1_pd(ntau);
    --v;
    a_buff -= ldr + 1;
    --work;

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp = _mm_setzero_pd();
        vd4_dtmp = _mm256_setzero_pd();

        /* v(1) = 1, start with the contribution of v(1) to tmp */
        vd2_dtmp3 = _mm_load1_pd(&a_buff[1 + j * ldr]);

        /* Compute tmp = c_A**T . v */
        for(k = 2; k <= (arows - 3); k += 4)
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);
        }
        if(k < arows)
        {
            /* load column elements of A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
        }
        /* Horizontal add of dtmp */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);
        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp3);

        /* Store the result in work */
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

        /* v(1) = 1, add the contribution of v(1) to c_A */
        vd2_inp = _mm_load_sd(&a_buff[1 + j * ldr]);
        vd2_inp = _mm_add_sd(vd2_inp, vd2_dtmp);
        _mm_storel_pd(&a_buff[1 + j * ldr], vd2_inp);

        for(k = 2; k <= (arows - 3); k += 4)
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
        }
        if(k < arows)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
            k += 2;
        }
        if(k == arows)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storel_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }
    }
}

/* v(m) is not stored explicitly and is assumed to be 1 */
__attribute__((aligned(512))) void fla_dlarf1l_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n,
                                                                     doublereal *a_buff,
                                                                     aocl_int64_t ldr,
                                                                     doublereal *v, doublereal ntau,
                                                                     doublereal *work)
{
    aocl_int64_t acols, arows, k_end;
    aocl_int64_t k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp, vd2_vj1, vd2_dtmp2;
    __m256d vd4_dtmp, vd4_inp, vd4_vj, vd4_dtmp2;
    __m128d vd2_ltmp, vd2_htmp;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    arows = m;
    acols = n;
    k_end = arows - 1;
    vd2_ntau = _mm_set1_pd(ntau);
    --v;
    a_buff -= ldr + 1;
    --work;

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp = _mm_setzero_pd();
        vd4_dtmp = _mm256_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (k_end - 3); k += 4)
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);
        }
        if(k < k_end)
        {
            /* load column elements of A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
            k += 2;
        }
        if(k == k_end)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* take dot product */
            vd2_dtmp2 = _mm_mul_pd(vd2_inp, vd2_vj1);
            vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp2);
        }
        /* Horizontal add of dtmp */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);
        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_htmp);
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        /* v(arows) = 1, add the contribution of v(arows) to tmp */
        vd2_dtmp = _mm_add_pd(vd2_dtmp, _mm_load1_pd(&a_buff[arows + j * ldr]));

        /* Store the result in work */
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
        for(k = 1; k <= (k_end - 3); k += 4)
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
        }
        if(k < k_end)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
            k += 2;
        }
        if(k == k_end)
        {
            /* load single remaining element from c_A and v */
            vd2_inp = _mm_load_sd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj1 = _mm_load_sd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd2_dtmp2 = _mm_mul_pd(vd2_dtmp, vd2_vj1);
            vd2_inp = _mm_add_pd(vd2_dtmp2, vd2_inp);
            _mm_storel_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }

        /* v(arows) = 1, add the contribution of v(arows) to c_A */
        vd2_inp = _mm_load_sd(&a_buff[arows + j * ldr]);
        vd2_inp = _mm_add_sd(vd2_inp, vd2_dtmp);
        _mm_storel_pd(&a_buff[arows + j * ldr], vd2_inp);
    }
}

#endif
