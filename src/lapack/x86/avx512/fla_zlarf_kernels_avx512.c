/******************************************************************************
 * * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

void fla_zlarf_left_apply_incv1_avx512(integer m, integer n, doublecomplex *a_buff, integer ldr,
                                       doublecomplex *v, doublecomplex *ntau, doublecomplex *work)
{
    integer acols, arows;
    integer k, j;
    __m128d vd2_inp;
    __m128d vd2_ntau, vd2_dtmp1, vd2_dtmp2, vd2_vj, vd2_vjr, vd2_vji;
    __m128d vd2_ltmp, vd2_htmp;
    __m256d vd4_inp, vd4_dtmp1, vd4_dtmp2, vd4_vj, vd4_vjr, vd4_vji;
    __m512d vd8_dtmp1, vd8_dtmp2, vd8_inp, vd8_vj, vd8_vjr, vd8_vji;
    __m256d vd4_ltmp, vd4_htmp;
    __m128d vd2_one_neg_one = _mm_set_pd(-1.0, 1.0);

    /* Apply the Householder rotation                      */
    /* on the rest of the matrix                           */
    /*    A = A - tau * v * v**T * A                       */
    /*      = A - v * tau * (A**T * v)**T                  */
    /* DGEMV and DGER operations are combined              */

    arows = m;
    acols = n;

    vd2_ntau = _mm_loadu_pd((const doublereal *)ntau);

    /* Compute A**T * v */
    for(j = 1; j <= acols; j++) /* for every column c_A of A */
    {
        vd2_dtmp1 = _mm_setzero_pd();
        vd2_dtmp2 = _mm_setzero_pd();
        vd4_dtmp1 = _mm256_setzero_pd();
        vd4_dtmp2 = _mm256_setzero_pd();
        vd8_dtmp1 = _mm512_setzero_pd();
        vd8_dtmp2 = _mm512_setzero_pd();

        /* Compute tmp = c_A**T . v */
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of A and v */
            /* Column is loaded in the following format */
            /* [Ar0, Ai0, Ar1, Ai1, Ar2, Ai2, Ar3, Ai3] */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            /* Vector v is loaded in following format */
            /* [Vr0, Vi0, Vr1, Vi1, Vr2, Vi2, Vr3, Vi3] */
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* Rearrange vd8 as follows */
            /* vd8_vjr = [Vr0, Vr0, Vr1, Vr1, Vr2, Vr2, Vr3, Vr3] */
            vd8_vjr = _mm512_permute_pd(vd8_vj, 0b00000000);
            /* vd8_vji = [Vi0, Vi0, Vi1, Vi1, Vi2, Vi2, Vi3, Vi3] */
            vd8_vji = _mm512_permute_pd(vd8_vj, 0b11111111);

            /* take dot product */
            vd8_dtmp1 = _mm512_fmadd_pd(vd8_inp, vd8_vjr, vd8_dtmp1);
            vd8_dtmp2 = _mm512_fmadd_pd(vd8_inp, vd8_vji, vd8_dtmp2);
        }
        if(k <= (arows - 1))
        {
            /* load column elements of A and v */
            /* Column loaded in the following format */
            /* [Ar0, Ai0, Ar1, Ai1] */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            /* Vector v is loaded in following format */
            /* [Vr0, Vi0, Vr1, Vi1] */
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            /* Rearrange vd4 as follows */
            /* vd4_vjr = [Vr0, Vr0, Vr1, Vr1] */
            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            /* vd4_vji = [Vi0, Vi0, Vi1, Vi1] */
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            /* take dot product */
            vd4_dtmp1 = _mm256_fmadd_pd(vd4_inp, vd4_vjr, vd4_dtmp1);
            vd4_dtmp2 = _mm256_fmadd_pd(vd4_inp, vd4_vji, vd4_dtmp2);

            k += 2;
        }
        if(k == arows)
        {
            /* load column elements of A and v */
            /* Column is loaded in following format */
            /* [Ar0, Ai0] */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);

            /* Vector c is loaded in following format */
            /* [Vr0, Vi0] */
            vd2_vj = _mm_loadu_pd((const doublereal *)&v[k]);

            /* Rearrange vd2 as follows */
            /* vd2vjr = [Vr0, Vr0] */
            vd2_vjr = _mm_permute_pd(vd2_vj, 0b00);
            /* vd2vji = [Vi0, Vi0] */
            vd2_vji = _mm_permute_pd(vd2_vj, 0b11);

            /* take dot product */
            vd2_dtmp1 = _mm_fmadd_pd(vd2_inp, vd2_vjr, vd2_dtmp1);
            vd2_dtmp2 = _mm_fmadd_pd(vd2_inp, vd2_vji, vd2_dtmp2);
            k += 1;
        }

        /* Reduce add the values in vd8_dtmp, vd4_dtmp and vd2_dtmp */

        /* Etract Upper and lower 256 bits of vd8_dtmp1 */
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp1);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp1, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp1 */
        vd4_dtmp1 = _mm256_add_pd(vd4_dtmp1, vd4_ltmp);
        vd4_dtmp1 = _mm256_add_pd(vd4_dtmp1, vd4_htmp);

        /* Etract Upper and lower 256 bits of vd8_dtmp2 */
        vd4_ltmp = _mm512_castpd512_pd256(vd8_dtmp2);
        vd4_htmp = _mm512_extractf64x4_pd(vd8_dtmp2, 0x1);

        /* Add the lower and upper 256 bits with vd4_dtmp2 */
        vd4_dtmp2 = _mm256_add_pd(vd4_dtmp2, vd4_ltmp);
        vd4_dtmp2 = _mm256_add_pd(vd4_dtmp2, vd4_htmp);

        /* Etract Upper and lower 128 bits of vd4_dtmp1 */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp1);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp1, 0x1);

        /* Add the lower and upper 128 bits and store in vd2_dtmp1 */
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_ltmp);
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_htmp);

        /* Etract Upper and lower 128 bits of vd4_dtmp2 */
        vd2_ltmp = _mm256_castpd256_pd128(vd4_dtmp2);
        vd2_htmp = _mm256_extractf128_pd(vd4_dtmp2, 0x1);

        /* Add the lower and upper 128 bits and store in vd2_dtmp2 */
        vd2_dtmp2 = _mm_add_pd(vd2_dtmp2, vd2_ltmp);
        vd2_dtmp2 = _mm_add_pd(vd2_dtmp2, vd2_htmp);

        /*
            Register vd2_dtmp1 = [ Ar * Vr, Ai * Vr ]
            Regsiter vd2_dtmp2 = [ Ar * Vi, Ai * Vi ]
            Since taking conjugate of A, the multiplcation result would
            bas as follows
            [ Ar * Vr + Ai * Vi, Ar * Vi - Ai * Vr ]
        */
        /* permuting vd2_dtmp2 = [ Ai * Vi, Ar * Vi ] */
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);
        /* multiplying vd2_dtmp2 with [1.0, -1.0] */
        /* vd2_dtmp1 = [ Ar * Vr, - Ai * Vr ] */
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        /* adding vd2_dtmp1 and vd2_dtmp2 */
        vd2_dtmp1 = _mm_add_pd(vd2_dtmp1, vd2_dtmp2);

        /* Store the result in work */
        _mm_storeu_pd((doublereal *)&work[j], vd2_dtmp1);

        /* Take conjugate of tmp */
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);

        /* Compute tmp = ntau * tmp */
        /* vd2_vjr = [ Tr, Tr ] */
        vd2_vjr = _mm_permute_pd(vd2_ntau, 0b00);
        /* vd2_vji = [ Ti. Ti ] */
        vd2_vji = _mm_permute_pd(vd2_ntau, 0b11);
        /*
            tmp = [ Kr, Ki ]
            tau = [ Tr, Ti ]
            Muliplication will be as follows
            [Kr * Tr - Ki * Ti, Kr * Ti + Ki * Tr]
        */
        /* vd2_temp = [ Kr, -Ki ] */
        vd2_dtmp2 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        /* Permute vd2_dtmp2 = [ -Ki, Kr ] */
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);
        /* vd2_dtmp1 = [ Kr * Tr, Ki * Tr ] */
        vd2_dtmp1 = _mm_mul_pd(vd2_dtmp1, vd2_vjr);
        /* vd2_dtmp1 = [Kr * Tr - Ki * Ti, Kr * Ti + Ki * Tr] */
        vd2_dtmp1 = _mm_fmadd_pd(vd2_vji, vd2_dtmp2, vd2_dtmp1);

        /* Set the first element as negative [ Pr -Pi ] */
        vd2_dtmp2 = _mm_mul_pd(vd2_dtmp1, vd2_one_neg_one);
        /* Shuffle vd2_dtmp2 = [ -Pi Pr ] */
        vd2_dtmp2 = _mm_permute_pd(vd2_dtmp2, 0b01);

        /* Broadcast value to __m512d and _m256d */
        vd4_dtmp1 = _mm256_castpd128_pd256(vd2_dtmp1);
        vd4_dtmp1 = _mm256_insertf128_pd(vd4_dtmp1, vd2_dtmp1, 0x1);

        vd4_dtmp2 = _mm256_castpd128_pd256(vd2_dtmp2);
        vd4_dtmp2 = _mm256_insertf128_pd(vd4_dtmp2, vd2_dtmp2, 0x1);

        vd8_dtmp1 = _mm512_castpd256_pd512(vd4_dtmp1);
        vd8_dtmp1 = _mm512_insertf64x4(vd8_dtmp1, vd4_dtmp1, 0x1);

        vd8_dtmp2 = _mm512_castpd256_pd512(vd4_dtmp2);
        vd8_dtmp2 = _mm512_insertf64x4(vd8_dtmp2, vd4_dtmp2, 0x1);

        /* Compute c_A + (tmp * v')' = c_A + (v * tmp')*/
        for(k = 1; k <= (arows - 3); k += 4)
        {
            /* load column elements of c_A and v */
            /* vd8_inp = [ Ar0, Ai0, Ar1, Ai1, Ar2, Ai2, Ar3, Ai3  ] */
            vd8_inp = _mm512_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            /* vd8_vj =  [ Vr0, Vi0, Vr1, Vi1, Vr2, Vi2, Vr3, Vi3 ] */
            vd8_vj = _mm512_loadu_pd((const doublereal *)&v[k]);

            /* vd8_vjr = [ Vr0, Vr0, Vr1, Vr1, Vr2, Vr2, Vr3, Vr3 ] */
            vd8_vjr = _mm512_permute_pd(vd8_vj, 0b00000000);
            /* vd8_vji = [ Vi0, Vi0, Vi1, Vi1, Vi2, Vi2, Vi3, Vi3 ] */
            vd8_vji = _mm512_permute_pd(vd8_vj, 0b11111111);

            /* mul by dtmp, add and store */
            /* inp + [  Pr * Vr, Pi * Vr ] */
            vd8_inp = _mm512_fmadd_pd(vd8_dtmp1, vd8_vjr, vd8_inp);
            /* inp + [ Pr * Vr - Pi * Vi, Pi * Vr + Pr * Vi ] */
            vd8_inp = _mm512_fmadd_pd(vd8_dtmp2, vd8_vji, vd8_inp);

            _mm512_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd8_inp);
        }
        if(k <= (arows - 1))
        {
            /* Same steps followed as in above loop */
            /* load column elements of c_A and v */
            vd4_inp = _mm256_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd4_vj = _mm256_loadu_pd((const doublereal *)&v[k]);

            vd4_vjr = _mm256_permute_pd(vd4_vj, 0b0000);
            vd4_vji = _mm256_permute_pd(vd4_vj, 0b1111);

            vd4_inp = _mm256_fmadd_pd(vd4_dtmp1, vd4_vjr, vd4_inp);
            vd4_inp = _mm256_fmadd_pd(vd4_dtmp2, vd4_vji, vd4_inp);

            _mm256_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd4_inp);
            k += 2;
        }
        if(k == arows)
        {
            /* load column elements of c_A and v */
            vd2_inp = _mm_loadu_pd((const doublereal *)&a_buff[k + j * ldr]);
            vd2_vj = _mm_loadu_pd((const doublereal *)&v[k]);

            vd2_vjr = _mm_permute_pd(vd2_vj, 0b00);
            vd2_vji = _mm_permute_pd(vd2_vj, 0b11);

            /* mul by dtmp, add and store */
            vd2_inp = _mm_fmadd_pd(vd2_dtmp1, vd2_vjr, vd2_inp);
            vd2_inp = _mm_fmadd_pd(vd2_dtmp2, vd2_vji, vd2_inp);

            _mm_storeu_pd((doublereal *)&a_buff[k + j * ldr], vd2_inp);
        }
    }
}
#endif