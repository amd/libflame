/******************************************************************************
 * * Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/**
 * Applies the Householder reflector for DLARF (incv == 1, plain variant).
 * Required: m >= 1
 */
void fla_dlarf_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                       doublereal *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp, vd2_ntau, vd2_ltmp, vd2_htmp;
    __m128d vd2_dtmp, vd2_dtmp2, vd2_vj1;
    __m256d vd4_inp, vd4_dtmp, vd4_vj, vd4_dtmp2;
    __m256d vd4_ltmp, vd4_htmp;
    __m512d vd8_dtmp, vd8_inp, vd8_vj, vd8_dtmp2;

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
        vd8_dtmp = _mm512_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 7); k += 8)
        {
            /* load column elements of A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd8_dtmp2 = _mm512_mul_pd(vd8_inp, vd8_vj);
            vd8_dtmp = _mm512_add_pd(vd8_dtmp, vd8_dtmp2);
        }
        if(k <= (arows - 3))
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);

            k += 4;
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

        /* Reduce add the values in vd8_dtmp, vd4_dtmp and vd2_dtmp*/

        /* Etract Upper and lower 256 bits of vd8_dtmp*/
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp */
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_ltmp);
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_htmp);

        /* Horizontal add of dtmp */
        vd4_dtmp = _mm256_hadd_pd(vd4_dtmp, vd4_dtmp);

        /* Etract Upper and lower 128 bits of vd4_dtmp*/
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);
        /* Add the lower and upper 128 bits and store in vd2_ltmp */
        vd2_ltmp = _mm_add_pd(vd2_htmp, vd2_ltmp);

        /* Horizontal add of vd2_ltmp and vd2_dtmp */
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);

        /* Store the result in work */
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        vd8_dtmp = _mm512_castpd256_pd512(vd4_dtmp);
        vd8_dtmp = _mm512_insertf64x4(vd8_dtmp, vd4_dtmp, 0x1);

        /* Compute c_A + tmp * v */
        for(k = 1; k <= (arows - 7); k += 8)
        {
            /* load column elements of c_A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd8_dtmp2 = _mm512_mul_pd(vd8_dtmp, vd8_vj);
            vd8_inp = _mm512_add_pd(vd8_dtmp2, vd8_inp);
            _mm512_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd8_inp);
        }
        if(k <= (arows - 3))
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
            k += 4;
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

/**
 * Applies the Householder reflector for DLARF1F (incv == 1, v(1) implicitly 1).
 * Required: m >= 1
 */
void fla_dlarf1f_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                         aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                         doublereal *work)
{
    aocl_int64_t acols, arows;
    aocl_int64_t k, j;
    __m128d vd2_inp, vd2_ntau, vd2_ltmp, vd2_htmp;
    __m128d vd2_dtmp, vd2_dtmp2, vd2_vj1, vd2_dtmp3;
    __m256d vd4_inp, vd4_dtmp, vd4_vj, vd4_dtmp2;
    __m256d vd4_ltmp, vd4_htmp;
    __m512d vd8_dtmp, vd8_inp, vd8_vj, vd8_dtmp2;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    /* v(1) is not stored explicitly and is assumed to be 1 */

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
        vd8_dtmp = _mm512_setzero_pd();

        /* v(1) = 1, start with the contribution of v(1) to tmp */
        vd2_dtmp3 = _mm_load1_pd(&a_buff[1 + j * ldr]);

        /* Compute tmp = c_A**T . v */
        for(k = 2; k <= (arows - 7); k += 8)
        {
            /* load column elements of A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd8_dtmp2 = _mm512_mul_pd(vd8_inp, vd8_vj);
            vd8_dtmp = _mm512_add_pd(vd8_dtmp, vd8_dtmp2);
        }
        if(k <= (arows - 3))
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);

            k += 4;
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

        /* Reduce add the values in vd8_dtmp, vd4_dtmp and vd2_dtmp*/

        /* Etract Upper and lower 256 bits of vd8_dtmp*/
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp */
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_ltmp);
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_htmp);

        /* Horizontal add of dtmp */
        vd4_dtmp = _mm256_hadd_pd(vd4_dtmp, vd4_dtmp);

        /* Etract Upper and lower 128 bits of vd4_dtmp*/
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);
        /* Add the lower and upper 128 bits and store in vd2_ltmp */
        vd2_ltmp = _mm_add_pd(vd2_htmp, vd2_ltmp);

        /* Horizontal add of vd2_ltmp and vd2_dtmp */
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_dtmp3);

        /* Store the result in work */
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        vd8_dtmp = _mm512_castpd256_pd512(vd4_dtmp);
        vd8_dtmp = _mm512_insertf64x4(vd8_dtmp, vd4_dtmp, 0x1);

        /* Compute c_A + tmp * v */

        /* v(1) = 1, add the contribution of v(1) to c_A */
        vd2_inp = _mm_load_sd(&a_buff[1 + j * ldr]);
        vd2_inp = _mm_add_sd(vd2_inp, vd2_dtmp);
        _mm_storel_pd(&a_buff[1 + j * ldr], vd2_inp);

        for(k = 2; k <= (arows - 7); k += 8)
        {
            /* load column elements of c_A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd8_dtmp2 = _mm512_mul_pd(vd8_dtmp, vd8_vj);
            vd8_inp = _mm512_add_pd(vd8_dtmp2, vd8_inp);
            _mm512_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd8_inp);
        }
        if(k <= (arows - 3))
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
            k += 4;
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

/**
 * Applies the Householder reflector for DLARF1L (incv == 1, v(m) implicitly 1).
 * Required: m >= 1
 */
void fla_dlarf1l_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                         aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                         doublereal *work)
{
    aocl_int64_t acols, arows, k_end;
    aocl_int64_t k, j;
    __m128d vd2_inp, vd2_ntau, vd2_ltmp, vd2_htmp;
    __m128d vd2_dtmp, vd2_dtmp2, vd2_vj1;
    __m256d vd4_inp, vd4_dtmp, vd4_vj, vd4_dtmp2;
    __m256d vd4_ltmp, vd4_htmp;
    __m512d vd8_dtmp, vd8_inp, vd8_vj, vd8_dtmp2;

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */
    /* v(m) is not stored explicitly and is assumed to be 1 */

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
        vd8_dtmp = _mm512_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (k_end - 7); k += 8)
        {
            /* load column elements of A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd8_dtmp2 = _mm512_mul_pd(vd8_inp, vd8_vj);
            vd8_dtmp = _mm512_add_pd(vd8_dtmp, vd8_dtmp2);
        }
        if(k <= (k_end - 3))
        {
            /* load column elements of A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* take dot product */
            vd4_dtmp2 = _mm256_mul_pd(vd4_inp, vd4_vj);
            vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_dtmp2);

            k += 4;
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

        /* Reduce add the values in vd8_dtmp, vd4_dtmp and vd2_dtmp*/

        /* Etract Upper and lower 256 bits of vd8_dtmp*/
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp */
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_ltmp);
        vd4_dtmp = _mm256_add_pd(vd4_dtmp, vd4_htmp);

        /* Horizontal add of dtmp */
        vd4_dtmp = _mm256_hadd_pd(vd4_dtmp, vd4_dtmp);

        /* Etract Upper and lower 128 bits of vd4_dtmp*/
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp, 0x1);
        /* Add the lower and upper 128 bits and store in vd2_ltmp */
        vd2_ltmp = _mm_add_pd(vd2_htmp, vd2_ltmp);

        /* Horizontal add of vd2_ltmp and vd2_dtmp */
        vd2_dtmp = _mm_hadd_pd(vd2_dtmp, vd2_dtmp);

        vd2_dtmp = _mm_add_pd(vd2_dtmp, vd2_ltmp);

        /* v(arows) = 1, add the contribution of v(arows) to tmp */
        vd2_dtmp = _mm_add_pd(vd2_dtmp, _mm_load1_pd(&a_buff[arows + j * ldr]));

        /* Store the result in work */
        _mm_storel_pd((doublereal *)&work[j], vd2_dtmp);

        /* Compute tmp = - tau * tmp */
        vd2_dtmp = _mm_mul_pd(vd2_dtmp, vd2_ntau);
        vd4_dtmp = _mm256_castpd128_pd256(vd2_dtmp);
        vd4_dtmp = _mm256_insertf128_pd(vd4_dtmp, vd2_dtmp, 0x1);

        vd8_dtmp = _mm512_castpd256_pd512(vd4_dtmp);
        vd8_dtmp = _mm512_insertf64x4(vd8_dtmp, vd4_dtmp, 0x1);

        /* Compute c_A + tmp * v */
        for(k = 1; k <= (k_end - 7); k += 8)
        {
            /* load column elements of c_A and v */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd8_dtmp2 = _mm512_mul_pd(vd8_dtmp, vd8_vj);
            vd8_inp = _mm512_add_pd(vd8_dtmp2, vd8_inp);
            _mm512_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd8_inp);
        }
        if(k <= (k_end - 3))
        {
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* mul by dtmp, add and store */
            vd4_dtmp2 = _mm256_mul_pd(vd4_dtmp, vd4_vj);
            vd4_inp = _mm256_add_pd(vd4_dtmp2, vd4_inp);
            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
            k += 4;
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