/* ../netlib/dlarf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/*
 *     Modifications Copyright (c) 2024-2025 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLA_f2c.h" /* Table of constant values */
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static aocl_int64_t c__1 = 1;
/* > \brief \b DLARF applies an elementary reflector to a general rectangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, LDC, M, N */
/* DOUBLE PRECISION TAU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARF applies a real elementary reflector H to a real m by n matrix */
/* > C, from either the left or the right. H is represented in the form */
/* > */
/* > H = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar and v is a real vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': form H * C */
/* > = 'R': form C * H */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension */
/* > (1 + (M-1)*f2c_dabs(INCV)) if SIDE = 'L' */
/* > or (1 + (N-1)*f2c_dabs(INCV)) if SIDE = 'R' */
/* > The vector v in the representation of H. V is not used if */
/* > TAU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* > INCV is INTEGER */
/* > The increment between elements of v. INCV <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (LDC,N) */
/* > On entry, the m by n matrix C. */
/* > On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* > or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension */
/* > (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dlarf_(char *side, aocl_int_t *m, aocl_int_t *n, doublereal *v, aocl_int_t *incv,
            doublereal *tau, doublereal *c__, aocl_int_t *ldc, doublereal *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlarf(side, m, n, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_dlarf(side, &m_64, &n_64, v, &incv_64, tau, c__, &ldc_64, work);
#endif
}

void aocl_lapack_dlarf(char *side, aocl_int64_t *m, aocl_int64_t *n, doublereal *v,
                       aocl_int64_t *incv, doublereal *tau, doublereal *c__, aocl_int64_t *ldc,
                       doublereal *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlarf inputs: side %c, m %" FLA_IS ", n %" FLA_IS ", incv %" FLA_IS
                      ", ldc %" FLA_IS "",
                      *side, *m, *n, *incv, *ldc);
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset;
    doublereal d__1;
    /* Local variables */
    aocl_int64_t i__;
    logical applyleft;
#ifdef FLA_ENABLE_AMD_OPT
    extern void fla_dlarf_small_incv1_simd(aocl_int64_t lastv, aocl_int64_t lastc, double *c__,
                                           aocl_int64_t ldc, double *v, double tau, double *work);
#if !FLA_ENABLE_AOCL_BLAS
#endif
    void fla_dlarf_left_tuning_params(aocl_int64_t m, aocl_int64_t n, FLA_Bool * use_blocked_flag,
                                      aocl_int64_t * nthreads);
    void fla_dlarf_right_tuning_params(aocl_int64_t m, aocl_int64_t n, aocl_int64_t * block_size,
                                       aocl_int64_t * nthreads);
#endif
#if !FLA_ENABLE_AOCL_BLAS
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
#endif
    aocl_int64_t lastc, lastv;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    applyleft = lsame_(side, "L", 1, 1);
    lastv = 0;
    lastc = 0;
    if(*tau != 0.)
    {
        /* Set up variables for scanning V. LASTV begins pointing to the end */
        /* of V. */
        if(applyleft)
        {
            lastv = *m;
        }
        else
        {
            lastv = *n;
        }
        if(*incv > 0)
        {
            i__ = (lastv - 1) * *incv + 1;
        }
        else
        {
            i__ = 1;
        }
        /* Look for the last non-zero row in V. */
        while(lastv > 0 && v[i__] == 0.)
        {
            --lastv;
            i__ -= *incv;
        }
        if(applyleft)
        {
            /* Scan for the last non-zero column in C(1:lastv,:). */
            lastc = aocl_lapack_iladlc(&lastv, n, &c__[c_offset], ldc);
        }
        else
        {
            /* Scan for the last non-zero row in C(:,1:lastv). */
            lastc = aocl_lapack_iladlr(m, &lastv, &c__[c_offset], ldc);
        }
    }
    /* Note that lastc.eq.0 renders the BLAS operations null;
    no special */
    /* case is needed at this level. */
    if(applyleft)
    {
        /* Form H * C */
        if(lastv > 0 && lastc > 0)
        {
            d__1 = -(*tau);

#ifndef FLA_ENABLE_AMD_OPT
            /* w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */
            aocl_blas_dgemv("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &v[1], incv,
                            &c_b5, &work[1], &c__1);

            /* C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T */
            aocl_blas_dger(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset],
                           ldc);
#else
            /* Get threshold sizes to take optimized path*/
            FLA_Bool min_lastc_lastv = (lastc <= FLA_DGEMV_DGER_SIMD_SMALL_THRESH)
                                       && (lastv >= FLA_DGEMV_DGER_SIMD_SMALL_THRESH_M
                                           && lastv <= FLA_DGEMV_DGER_SIMD_SMALL_THRESH);

            /* Initialize global context data */
            aocl_fla_init();

            /* If the size of the matrix is small and incv =1, use the optimized path */
            if(min_lastc_lastv && *incv == c__1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
            {
                /* Call optimized routine */
                fla_dlarf_small_incv1_simd(lastv, lastc, c__, *ldc, v, d__1, work);
            }
            else
            {

                FLA_Bool use_blocked = 0;
                aocl_int64_t opt_nthreads = 1;

                fla_dlarf_left_tuning_params(lastv, lastc, &use_blocked, &opt_nthreads);

                /* If use_blocked is 1, process in blocks */
                if(use_blocked)
                {
                    /* Process in blocks */
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel for num_threads(opt_nthreads) private(i__)
#endif
                    /* Loop for each column of C */
                    for(i__ = 1; i__ <= lastc; ++i__)
                    {
                        /* W(i) =  C(1:lastv,i)**T * v(1:lastv,1)  */
                        work[i__]
                            = aocl_blas_ddot(&lastv, &v[1], incv, &c__[i__ * *ldc + 1], &c__1);
                        /* C(1:lastv,i) = C(1:lastv,i) - v(1:lastv,1) * -tau * W(i) */
                        doublereal d__2 = -(*tau) * work[i__];
                        aocl_blas_daxpy(&lastv, &d__2, &v[1], incv, &c__[i__ * *ldc + 1], &c__1);
                    }
                }
                else
                {
                    /* Process in a single call */
                    /* w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */
#if FLA_ENABLE_AOCL_BLAS && defined(BLIS_KERNELS_ZEN4)
                    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *incv > 0)
                    {
                        /* Use direct single threaded BLIS kernel */
                        bli_dgemv_t_zen4_int(BLIS_CONJUGATE, BLIS_NO_CONJUGATE, lastv, lastc, &c_b4,
                                             &c__[c_offset], 1, *ldc, &v[1], *incv, &c_b5, &work[1],
                                             c__1, NULL);
                    }
                    else
                    {
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp teams num_teams(1) thread_limit(1)
#endif
                        {
                            aocl_blas_dgemv("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc,
                                            &v[1], incv, &c_b5, &work[1], &c__1);
                        }
                    }
#else
                    aocl_blas_dgemv("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &v[1],
                                    incv, &c_b5, &work[1], &c__1);
#endif
                    /* C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T*/
                    aocl_blas_dger(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1,
                                   &c__[c_offset], ldc);
                }
            }
#endif /* FLA_ENABLE_AMD_OPT */
        }
    }
    else
    {
        /* Form C * H */
        if(lastv > 0)
        {
#ifndef FLA_ENABLE_AMD_OPT
            /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
            aocl_blas_dgemv("No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc, &v[1], incv,
                            &c_b5, &work[1], &c__1);
            /* C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T */
            d__1 = -(*tau);
            aocl_blas_dger(&lastc, &lastv, &d__1, &work[1], &c__1, &v[1], incv, &c__[c_offset],
                           ldc);
#else
            aocl_int64_t opt_nthreads = 1;
            aocl_int64_t nb = 0;

            fla_dlarf_right_tuning_params(lastc, lastv, &nb, &opt_nthreads);

            d__1 = -(*tau);

            /* If nb is non zero, process in blocks */
            if(nb && opt_nthreads > 1)
            {
                /* The first panel will process starting unaligned elements
                 * to ensure that all other panels aligned memory addresses
                 */
                uint64_t unaligned_bytes
                    = ((uint64_t)(c__ + c_offset)) % ((nb * sizeof(doublereal)));
                aocl_int64_t first_thread_rows
                    = fla_min((nb - (unaligned_bytes / sizeof(doublereal))), lastc);
                aocl_int64_t panels
                    = (!!first_thread_rows) + (((lastc - first_thread_rows) + (nb - 1)) / nb);

#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel for num_threads(opt_nthreads) private(i__)
#endif
                for(i__ = 1; i__ <= panels; i__ += 1)
                {
                    aocl_int64_t completed_rows = i__ == 1 ? 0 : (i__ - 2) * nb + first_thread_rows;
                    aocl_int64_t current_block_size
                        = i__ == 1 ? first_thread_rows : fla_min(nb, lastc - completed_rows);
                    aocl_int64_t cur_idx = completed_rows + 1;
                    /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
                    aocl_blas_dgemv("No transpose", &current_block_size, &lastv, &c_b4,
                                    &c__[c_dim1 + cur_idx], ldc, &v[1], incv, &c_b5, &work[cur_idx],
                                    &c__1);
                    aocl_blas_dger(&current_block_size, &lastv, &d__1, &work[cur_idx], &c__1, &v[1],
                                   incv, &c__[c_dim1 + cur_idx], ldc);
                }
            }
            else
            {
                /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
#if FLA_ENABLE_AOCL_BLAS && defined(BLIS_KERNELS_ZEN4)
                aocl_fla_init();
                if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *incv > 0)
                {
                    bli_dgemv_n_zen4_int_40x2_st(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, lastc, lastv,
                                                 &c_b4, &c__[c_offset], c__1, c_dim1, &v[1], *incv,
                                                 &c_b5, &work[1], c__1, NULL);
                }
                else
                {
                    aocl_blas_dgemv("No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc,
                                    &v[1], incv, &c_b5, &work[1], &c__1);
                }
#else
                aocl_blas_dgemv("No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc, &v[1],
                                incv, &c_b5, &work[1], &c__1);
#endif
                /* C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T */
                d__1 = -(*tau);
                aocl_blas_dger(&lastc, &lastv, &d__1, &work[1], &c__1, &v[1], incv, &c__[c_offset],
                               ldc);
            }
#endif /* FLA_ENABLE_AMD_OPT */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLARF */
}
/* dlarf_ */

#ifdef FLA_ENABLE_AMD_OPT
void fla_dlarf_left_tuning_params(aocl_int64_t m, aocl_int64_t n, FLA_Bool *use_blocked_flag,
                                  aocl_int64_t *nthreads)
{
    extern int fla_thread_get_num_threads(void);
    aocl_int64_t num_elems = m * n;
    if(num_elems < FLA_DLARF_L_THRESH_UNBLOCKED)
    {
        *use_blocked_flag = 0;
        return;
    }

    aocl_int64_t max_available_threads = fla_thread_get_num_threads();

    /* Special case for 1 thread */
    if(max_available_threads == 1)
    {
        *nthreads = 1;
        /* Use blocked code for large sizes as it is more cache friendly */
        if(m > FLA_DLARF_L_ST_BLOCKED_THRESH_M && n > FLA_DLARF_L_ST_BLOCKED_THRESH_N)
        {
            *use_blocked_flag = 1;
        }
        else
        {
            *use_blocked_flag = 0;
        }
        return;
    }

    /* General case */

    aocl_int64_t opt_n_threads = 1;

    if(num_elems < FLA_DLARF_L_THRESH_THREAD_8)
    {
        opt_n_threads = fla_min(8, n / 2);
    }
    else if(num_elems < FLA_DLARF_L_THRESH_THREAD_64)
    {
        opt_n_threads = fla_min(64, n / 2);
    }
    else
    {
        opt_n_threads = fla_min(128, n / 2);
    }

    *use_blocked_flag = 1;
    *nthreads = fla_min(opt_n_threads, max_available_threads);
}

void fla_dlarf_right_tuning_params(aocl_int64_t m, aocl_int64_t n, aocl_int64_t *block_size,
                                   aocl_int64_t *nthreads)
{
    extern int fla_thread_get_num_threads(void);
    aocl_int64_t num_elems = m * n;

    if(num_elems < FLA_DLARF_R_THRESH_UNBLOCKED)
    {
        *block_size = 0;
        return;
    }

    aocl_int64_t max_available_threads = fla_thread_get_num_threads();
    aocl_int64_t opt_n_threads = 1;

    if(num_elems < FLA_DLARF_R_THRESH_THREAD_8)
    {
        opt_n_threads = fla_min(8, m / 2);
    }
    else if(num_elems < FLA_DLARF_R_THRESH_THREAD_64)
    {
        opt_n_threads = fla_min(64, m / 2);
    }
    else
    {
        opt_n_threads = fla_min(128, m / 2);
    }

    *block_size = FLA_DLARF_R_BLOCK_SIZE;
    *nthreads = fla_min(opt_n_threads, max_available_threads);
}

#endif /* FLA_ENABLE_AMD_OPT */
