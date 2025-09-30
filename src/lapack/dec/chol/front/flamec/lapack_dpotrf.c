/*
 * Copyright (c) 2021-2025 Advanced Micro Devices, Inc. All rights reserved.
 */

/* dpotrf.f -- translated by f2c and slightly modified */

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#if FLA_ENABLE_AMD_OPT
int fla_dpotrf_small_avx2(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, aocl_int64_t *info);
#endif

/* Threshold values used for tuning thread binding and workload partitioning.*/
#define PROC_BIND_CLOSE_THREADS 96
#define PROC_BIND_CLOSE_SIZE 9760

extern int fla_thread_get_num_threads();

void dpotrf_auto_tune_params(aocl_int64_t n, int *num_threads, aocl_int64_t *block_size)
{
    *block_size = FLA_POTRF_BLOCK_SIZE; // Default block size

    // Get maximum available threads
    int max_threads = fla_thread_get_num_threads();

    if(n <= 300)
    {
        *num_threads = 8;
        *block_size = 64;
    }
    else if(n <= 1024)
    {
        *num_threads = 8;
        *block_size = 128;
    }
    else if(n <= 1280)
    {
        *num_threads = 32;
        *block_size = 128;
    }
    else if(n <= 2048)
    {
        *num_threads = 32;
        *block_size = 128;
    }
    else if(n >= 9760)
    {
        // Large problems
        *num_threads = 128;
        *block_size = 256;
    }
    else
    {
        // Medium-large problems
        *num_threads = 64;
        *block_size = 224;
    }
    // Ensure we don't exceed available threads
    *num_threads = fla_min(*num_threads, max_threads);
}

/* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static doublereal c_b13 = -1.;
static doublereal c_b14 = 1.;

/* Subroutine */ int lapack_dpotrf(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                                   aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    aocl_int64_t j, jb, nb;
    logical upper;

#ifndef FLA_ENABLE_AOCL_BLAS
    logical lsame_(char *ca, char *cb, aocl_int64_t a, aocl_int64_t b);
#endif
    int lapack_dpotf2(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, aocl_int64_t *info);

    /*  DPOTRF computes the Cholesky factorization of a real symmetric */
    /*  positive definite matrix A. */

    /*  The factorization has the form */
    /*     A = U**T * U,  if UPLO = 'U', or */
    /*     A = L  * L**T,  if UPLO = 'L', */
    /*  where U is an upper triangular matrix and L is lower triangular. */

    /*  This is the block version of the algorithm, calling Level 3 BLAS. */

    /*  Arguments */
    /*  ========= */

    /*  UPLO    (input) CHARACTER*1 */
    /*          = 'U':  Upper triangle of A is stored; */
    /*          = 'L':  Lower triangle of A is stored. */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A.  N >= 0. */

    /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
    /*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
    /*          N-by-N upper triangular part of A contains the upper */
    /*          triangular part of the matrix A, and the strictly lower */
    /*          triangular part of A is not referenced.  If UPLO = 'L', the */
    /*          leading N-by-N lower triangular part of A contains the lower */
    /*          triangular part of the matrix A, and the strictly upper */
    /*          triangular part of A is not referenced. */

    /*          On exit, if INFO = 0, the factor U or L from the Cholesky */
    /*          factorization A = U**T*U or A = L*L**T. */

    /*  LDA     (input) INTEGER */
    /*          The leading dimension of the array A.  LDA >= fla_max(1,N). */

    /*  INFO    (output) INTEGER */
    /*          = 0:  successful exit */
    /*          < 0:  if INFO = -i, the i-th argument had an illegal value */
    /*          > 0:  if INFO = i, the leading minor of order i is not */
    /*                positive definite, and the factorization could not be */
    /*                completed. */

    /*  ===================================================================== */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
#ifndef FLA_ENABLE_AMD_OPT
#if AOCL_FLA_PROGRESS_H
    AOCL_FLA_PROGRESS_VAR;
    progress_step_count = 0;
#ifndef FLA_ENABLE_WINDOWS_BUILD
    if(!aocl_fla_progress_ptr)
        aocl_fla_progress_ptr = aocl_fla_progress;
#endif
#endif
#endif
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DPOTRF", &i__1, (ftnlen)6);
        return 0;
    }

    /*     Quick return if possible */

    if(*n == 0)
    {
        return 0;
    }

    /*     Determine the block size for this environment. */

    nb = aocl_lapack_ilaenv(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1);
    if(nb <= 1 || nb >= *n)
    {

        /*        Use unblocked code. */
#if FLA_ENABLE_AMD_OPT
        fla_dpotrf_small_avx2(uplo, n, &a[a_offset], lda, info);
#else
        lapack_dpotf2(uplo, n, &a[a_offset], lda, info);
#endif
    }
    else
    {
        /*        Use blocked code. */
        if(upper)
        {

            /*           Compute the Cholesky factorization A = U'*U. */

            i__1 = *n;
            i__2 = nb;

            for(j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {

                /*              Update and factorize the current diagonal block and test */
                /*              for non-positive-definiteness. */

                /* Computing MIN */
                i__3 = nb, i__4 = *n - j + 1;
                jb = fla_min(i__3, i__4);
                i__3 = j - 1;
#ifndef FLA_ENABLE_AMD_OPT
#if AOCL_FLA_PROGRESS_H
                if(aocl_fla_progress_ptr)
                {
                    progress_step_count += jb;
                    AOCL_FLA_PROGRESS_FUNC_PTR("DPOTRF", 6, &progress_step_count,
                                               &progress_thread_id, &progress_total_threads);
                }
#endif
#endif
                aocl_blas_dsyrk("Upper", "Transpose", &jb, &i__3, &c_b13, &a[j * a_dim1 + 1], lda, &c_b14,
                       &a[j + j * a_dim1], lda);
#if FLA_ENABLE_AMD_OPT
                fla_dpotrf_small_avx2("Upper", &jb, &a[j + j * a_dim1], lda, info);
#else
                lapack_dpotf2("Upper", &jb, &a[j + j * a_dim1], lda, info);
#endif
                if(*info != 0)
                {
                    goto L30;
                }
                if(j + jb <= *n)
                {

                    /* Compute the current block row. */

                    i__3 = *n - j - jb + 1;
                    i__4 = j - 1;

                    aocl_blas_dgemm("Transpose", "No transpose", &jb, &i__3, &i__4, &c_b13,
                           &a[j * a_dim1 + 1], lda, &a[(j + jb) * a_dim1 + 1], lda, &c_b14,
                           &a[j + (j + jb) * a_dim1], lda);
                    i__3 = *n - j - jb + 1;

                    aocl_blas_dtrsm("Left", "Upper", "Transpose", "Non-unit", &jb, &i__3, &c_b14,
                           &a[j + j * a_dim1], lda, &a[j + (j + jb) * a_dim1], lda);
                }
                /* L10: */
            }
        }
        else
        {

            /*           Compute the Cholesky factorization A = L*L'. */

            i__2 = *n;
            i__1 = nb;

            for(j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
            {

                /*              Update and factorize the current diagonal block and test */
                /*              for non-positive-definiteness. */

                /* Computing MIN */
                i__3 = nb, i__4 = *n - j + 1;
                jb = fla_min(i__3, i__4);
                i__3 = j - 1;
#ifndef FLA_ENABLE_AMD_OPT
#if AOCL_FLA_PROGRESS_H
                if(aocl_fla_progress_ptr)
                {
                    progress_step_count += jb;
                    AOCL_FLA_PROGRESS_FUNC_PTR("DPOTRF", 6, &progress_step_count,
                                               &progress_thread_id, &progress_total_threads);
                }
#endif
#endif
                aocl_blas_dsyrk("Lower", "No transpose", &jb, &i__3, &c_b13, &a[j + a_dim1], lda, &c_b14,
                       &a[j + j * a_dim1], lda);
#if FLA_ENABLE_AMD_OPT
                fla_dpotrf_small_avx2("Lower", &jb, &a[j + j * a_dim1], lda, info);
#else
                lapack_dpotf2("Lower", &jb, &a[j + j * a_dim1], lda, info);
#endif
                if(*info != 0)
                {
                    goto L30;
                }
                if(j + jb <= *n)
                {

                    /*                 Compute the current block column. */

                    i__3 = *n - j - jb + 1;
                    i__4 = j - 1;
                    aocl_blas_dgemm("No transpose", "Transpose", &i__3, &jb, &i__4, &c_b13,
                           &a[j + jb + a_dim1], lda, &a[j + a_dim1], lda, &c_b14,
                           &a[j + jb + j * a_dim1], lda);
                    i__3 = *n - j - jb + 1;
                    aocl_blas_dtrsm("Right", "Lower", "Transpose", "Non-unit", &i__3, &jb, &c_b14,
                           &a[j + j * a_dim1], lda, &a[j + jb + j * a_dim1], lda);
                }
                /* L20: */
            }
        }
    }
    goto L40;

L30:
    *info = *info + j - 1;

L40:
    return 0;

    /*     End of DPOTRF */

} /* dpotrf_ */

#ifdef FLA_OPENMP_MULTITHREADING

static size_t A21;
static size_t A12;
static size_t A22;

#define A(m, n, gm, gn, mb, nb, uplo) \
    (double *)get_tile_addr_triangle(A, gm, m, n, gm, gn, mb, nb, uplo)

static inline aocl_int64_t get_tile_lda(aocl_int64_t mb, aocl_int64_t gm, aocl_int64_t k)
{
    if(k < gm / mb)
        return mb;
    else
        return gm % mb;
}

static inline void *get_tile_addr_triangle(double *A, aocl_int64_t N, aocl_int64_t m, aocl_int64_t n, aocl_int64_t gm,
                                           aocl_int64_t gn, aocl_int64_t mb, aocl_int64_t nb, char *uplo)
{
    size_t eltsize = sizeof(double);
    size_t offset = 0;
    aocl_int64_t lm1 = N / mb;
    aocl_int64_t ln1 = N / nb;
    if(m < lm1)
    {
        if(n < ln1)
        {
            if(*uplo == 'U')
            {
                offset = mb * nb * (m + (n * (n + 1)) / 2);
            }
            else
            {
                offset = mb * nb * ((m - n) + (n * (2 * lm1 - (n - 1))) / 2);
            }
        }
        else
        {
            offset = A12 + ((size_t)mb * (gn % nb) * m);
        }
    }
    else
    {
        if(n < ln1)
        {
            offset = A21 + ((size_t)nb * (gm % mb) * n);
        }
        else
        {
            offset = A22;
        }
    }
    return (void *)((char *)A + (offset * eltsize));
}

static inline aocl_int64_t get_tile_rows(aocl_int64_t mt, aocl_int64_t mb, aocl_int64_t m, aocl_int64_t k)
{
    if(k < mt - 1)
        return mb;
    else if((m) % mb == 0)
        return mb;
    else
        return (m) % mb;
}

void dpotrf_tile(char *uplo, aocl_int64_t n, double *A, aocl_int64_t lda, aocl_int64_t *iinfo
#if AOCL_FLA_PROGRESS_H
                 ,
                 aocl_fla_progress_callback aocl_fla_progress_ptr, aocl_int64_t progress_step_count,
                 aocl_int64_t progress_total_threads
#endif
)
{
#pragma omp task depend(inout : A[0 : lda * n])
    {
        __builtin_prefetch(A, 1, 3);
#if AOCL_FLA_PROGRESS_H
        aocl_int64_t thread_id = omp_get_thread_num();
        if(aocl_fla_progress_ptr)
        {
            AOCL_FLA_PROGRESS_FUNC_PTR("DPOTRF", 6, &progress_step_count, &thread_id,
                                       &progress_total_threads);
        }
#endif
        lapack_dpotrf(uplo, &n, A, &lda, iinfo);
    }
}

void dtrsm_tile(char *side, char *uplo, char *transa, char *diag, aocl_int64_t m, aocl_int64_t n,
                double *alpha, double *A, aocl_int64_t lda, double *B, aocl_int64_t ldb)
{
    aocl_int64_t ak;
    if(*side == 'L')
        ak = m;
    else
        ak = n;
#pragma omp task depend(in : A[0 : (lda) * ak]) depend(inout : B[0 : (ldb) * (n)])
    {
        __builtin_prefetch(B, 1, 3);
        aocl_blas_dtrsm(side, uplo, transa, diag, &m, &n, alpha, A, &lda, B, &ldb);
    }
}

void dsyrk_tile(char *uplo, char *trans, aocl_int64_t n, aocl_int64_t k, double *alpha, double *A,
                aocl_int64_t lda, double *beta, double *C, aocl_int64_t ldc)
{
    aocl_int64_t ak;
    if(*trans == 'N')
        ak = k;
    else
        ak = n;
#pragma omp task depend(in : A[0 : (lda) * ak]) depend(inout : C[0 : (ldc) * (n)])
    {
        __builtin_prefetch(A, 1, 3);
        aocl_blas_dsyrk(uplo, trans, &n, &k, alpha, A, &lda, beta, C, &ldc);
    }
}

void dgemm_tile(char *transa, char *transb, aocl_int64_t m, aocl_int64_t n, aocl_int64_t k, double *alpha,
                double *A, aocl_int64_t lda, double *B, aocl_int64_t ldb, double *beta, double *C,
                aocl_int64_t ldc)
{
    aocl_int64_t ak;
    if(*transa == 'N')
        ak = k;
    else
        ak = m;

    aocl_int64_t bk;
    if(*transb == 'N')
        bk = n;
    else
        bk = k;
#pragma omp task depend(in : A[0 : (lda) * ak]) depend(in : B[0 : (ldb) * bk]) \
    depend(inout : C[0 : (ldc) * (n)])
    {
        __builtin_prefetch(C, 1, 3);
        aocl_blas_dgemm(transa, transb, &m, &n, &k, alpha, A, &lda, B, &ldb, beta, C, &ldc);
    }
}

void omp_dpotrf(char *uplo, double *A, aocl_int64_t *n, aocl_int64_t *lda, aocl_int64_t mt, aocl_int64_t mb, aocl_int64_t nb,
                aocl_int64_t gm, aocl_int64_t gn, aocl_int64_t *info)
{
#if AOCL_FLA_PROGRESS_H
    AOCL_FLA_PROGRESS_VAR;
    progress_step_count = 0;
#ifndef FLA_ENABLE_WINDOWS_BUILD
    if(!aocl_fla_progress_ptr)
        aocl_fla_progress_ptr = aocl_fla_progress;
#endif
#endif
    if(*uplo == 'L')
    {
        /*
         * Lower triangular factorization: A = L * L^T
         *
         * Visual representation of tile dependencies:
         *
         * k=0: [L00]
         *      [L10] <- dtrsm(L00, L10)
         *      [L20] <- dtrsm(L00, L20)
         *      [L30] <- dtrsm(L00, L30)
         *
         * k=1: [L00] [   ]
         *      [L10] [L11] <- dpotrf(L11 - L10*L10^T)
         *      [L20] [L21] <- dtrsm(L11, L21)
         *      [L30] [L31] <- dtrsm(L11, L31)
         *
         * k=2: [L00] [   ] [   ]
         *      [L10] [L11] [   ]
         *      [L20] [L21] [L22] <- dpotrf(L22 - L20*L20^T - L21*L21^T)
         *      [L30] [L31] [L32] <- dtrsm(L22, L32)
         */

        for(aocl_int64_t k = 0; k < mt; k++)
        {
            aocl_int64_t mvak = get_tile_rows(mt, mb, *n, k);
            aocl_int64_t ldak = get_tile_lda(mb, gm, k);
#if AOCL_FLA_PROGRESS_H
            progress_step_count += mvak;
#endif
            // Factorize diagonal tile: A(k,k) = L(k,k) * L(k,k)^T
            dpotrf_tile(uplo, mvak, A(k, k, gm, gn, mb, nb, uplo), ldak, info
#if AOCL_FLA_PROGRESS_H
                        ,
                        aocl_fla_progress_ptr, progress_step_count, progress_total_threads
#endif
            );
            // Update column tiles below diagonal: A(m,k) = L(m,k) where m > k
            for(aocl_int64_t m = k + 1; m < mt; m++)
            {
                aocl_int64_t mvam = get_tile_rows(mt, mb, *n, m);
                aocl_int64_t ldam = get_tile_lda(mb, gm, m);
                dtrsm_tile("R", "L", "T", "N", mvam, mb, &c_b14, A(k, k, gm, gn, mb, nb, uplo),
                           ldak, A(m, k, gm, gn, mb, nb, uplo), ldam);
            }

            // Update trailing submatrix: A(m,n) -= L(m,k) * L(n,k)^T
            for(aocl_int64_t m = k + 1; m < mt; m++)
            {
                aocl_int64_t mvam = get_tile_rows(mt, mb, *n, m);
                aocl_int64_t ldam = get_tile_lda(mb, gm, m);

                // Update diagonal tile: A(m,m) -= L(m,k) * L(m,k)^T
                dsyrk_tile("L", "N", mvam, mb, &c_b13, A(m, k, gm, gn, mb, nb, uplo), ldam, &c_b14,
                           A(m, m, gm, gn, mb, nb, uplo), ldam);

                // Update off-diagonal tiles: A(m,n) -= L(m,k) * L(n,k)^T
                for(aocl_int64_t n = k + 1; n < m; n++)
                {
                    aocl_int64_t ldan = get_tile_lda(mb, gm, n);
                    dgemm_tile("N", "Transpose", mvam, mb, mb, &c_b13,
                               A(m, k, gm, gn, mb, nb, uplo), ldam, A(n, k, gm, gn, mb, nb, uplo),
                               ldan, &c_b14, A(m, n, gm, gn, mb, nb, uplo), ldam);
                }
            }
        }
    }
    else
    {
        /*
         * Upper triangular factorization: A = U^T * U
         *
         * Visual representation of tile dependencies:
         *
         * k=0: [U00] [U01] [U02] [U03]
         *      [   ] [   ] [   ] [   ]
         *      [   ] [   ] [   ] [   ]
         *      [   ] [   ] [   ] [   ]
         *      After k=0: U01, U02, U03 solved via dtrsm
         *
         * k=1: [U00] [U01] [U02] [U03]
         *      [   ] [U11] [U12] [U13]
         *      [   ] [   ] [   ] [   ]
         *      [   ] [   ] [   ] [   ]
         *      U11 = dpotrf(A11 - U01^T*U01), U12,U13 solved via dtrsm
         *
         * k=2: [U00] [U01] [U02] [U03]
         *      [   ] [U11] [U12] [U13]
         *      [   ] [   ] [U22] [U23]
         *      [   ] [   ] [   ] [   ]
         *      U22 = dpotrf(A22 - U02^T*U02 - U12^T*U12), U23 solved via dtrsm
         */

        for(aocl_int64_t k = 0; k < mt; k++)
        {
            aocl_int64_t nvak = get_tile_rows(mt, mb, *n, k);
            aocl_int64_t ldak = get_tile_lda(mb, gm, k);

            // Factorize diagonal tile: A(k,k) = U(k,k)^T * U(k,k)
            dpotrf_tile(uplo, nvak, A(k, k, gm, gn, mb, nb, uplo), ldak, info
#if AOCL_FLA_PROGRESS_H
                        ,
                        aocl_fla_progress_ptr, progress_step_count, progress_total_threads
#endif
            );
#if AOCL_FLA_PROGRESS_H
            if(aocl_fla_progress_ptr)
            {
                progress_step_count += nvak;
                AOCL_FLA_PROGRESS_FUNC_PTR("DPOTRF", 6, &progress_step_count, &progress_thread_id,
                                           &progress_total_threads);
            }
#endif
            // Update row tiles to the right of diagonal: A(k,m) = U(k,m) where m > k
            for(aocl_int64_t m = k + 1; m < mt; m++)
            {
                aocl_int64_t nvam = get_tile_rows(mt, mb, *n, m);
                dtrsm_tile("L", "U", "T", "N", mb, nvam, &c_b14, A(k, k, gm, gn, mb, nb, uplo),
                           ldak, A(k, m, gm, gn, mb, nb, uplo), ldak);
            }

            // Update trailing submatrix: A(n,m) -= U(k,n)^T * U(k,m)
            for(aocl_int64_t m = k + 1; m < mt; m++)
            {
                aocl_int64_t nvam = get_tile_rows(mt, mb, *n, m);
                aocl_int64_t ldam = get_tile_lda(mb, gm, m);

                // Update diagonal tile: A(m,m) -= U(k,m)^T * U(k,m)
                dsyrk_tile("U", "T", nvam, mb, &c_b13, A(k, m, gm, gn, mb, nb, uplo), ldak, &c_b14,
                           A(m, m, gm, gn, mb, nb, uplo), ldam);

                // Update off-diagonal tiles: A(n,m) -= U(k,n)^T * U(k,m)
                for(aocl_int64_t n = k + 1; n < m; n++)
                {
                    aocl_int64_t ldan = get_tile_lda(mb, gm, n);
                    dgemm_tile("T", "N", mb, nvam, mb, &c_b13, A(k, n, gm, gn, mb, nb, uplo), ldak,
                               A(k, m, gm, gn, mb, nb, uplo), ldak, &c_b14,
                               A(n, m, gm, gn, mb, nb, uplo), ldan);
                }
            }
        }
    }
}

void dlacpy_tile(aocl_int64_t m, aocl_int64_t n, double *A, aocl_int64_t lda, double *B, aocl_int64_t ldb)
{
#pragma omp task depend(in : A[0 : (lda) * (n)]) depend(out : B[0 : (ldb) * (n)])
    {
        aocl_lapack_dlacpy("Full", &m, &n, A, &lda, B, &ldb);
    }
}

/* Convert column-major (CM) matrix layout to tiled (CCRB).*/
void matrix_tile(double *pA, aocl_int64_t lda, double *A, aocl_int64_t nb, aocl_int64_t mb, char *uplo, aocl_int64_t mt,
                 aocl_int64_t nt, aocl_int64_t gm, aocl_int64_t gn, aocl_int64_t N)
{
    for(aocl_int64_t m = 0; m < mt; m++)
    {
        aocl_int64_t ldt = get_tile_lda(mb, gm, m);
        aocl_int64_t n_start = (*uplo == 'U' ? m : 0);
        aocl_int64_t n_end = (*uplo == 'U' ? nt : m + 1);
        for(aocl_int64_t n = n_start; n < n_end; n++)
        {
            aocl_int64_t x1 = 0;
            aocl_int64_t y1 = 0;
            aocl_int64_t x2 = n == nt - 1 ? (N - 1) % nb + 1 : nb;
            aocl_int64_t y2 = m == mt - 1 ? (N - 1) % mb + 1 : mb;

            double *f77 = &pA[(size_t)nb * lda * n + (size_t)mb * m];
            double *bdl = (double *)get_tile_addr_triangle(A, N, m, n, gm, gn, nb, nb, uplo);
            dlacpy_tile(y2 - y1, x2 - x1, &(f77[x1 * lda + y1]), lda, &(bdl[x1 * nb + y1]), ldt);
        }
    }
}
/* Convert tiled (CCRB) to column-major (CM) matrix layout.
   Out-of-place. */
void matrix_untile(double *pA, aocl_int64_t lda, double *A, aocl_int64_t nb, aocl_int64_t mb, char *uplo,
                   aocl_int64_t mt, aocl_int64_t nt, aocl_int64_t gm, aocl_int64_t gn, aocl_int64_t N)
{
    for(aocl_int64_t m = 0; m < mt; m++)
    {
        aocl_int64_t ldt = get_tile_lda(mb, gm, m);
        aocl_int64_t n_start = (*uplo == 'U' ? m : 0);
        aocl_int64_t n_end = (*uplo == 'U' ? nt : m + 1);
        for(aocl_int64_t n = n_start; n < n_end; n++)
        {
            aocl_int64_t x1 = 0;
            aocl_int64_t y1 = 0;
            aocl_int64_t x2 = n == nt - 1 ? (N - 1) % nb + 1 : nb;
            aocl_int64_t y2 = m == mt - 1 ? (N - 1) % mb + 1 : mb;
            double *f77 = &pA[(size_t)nb * lda * n + (size_t)mb * m];
            double *bdl = (double *)get_tile_addr_triangle(A, N, m, n, gm, gn, nb, nb, uplo);

            dlacpy_tile(y2 - y1, x2 - x1, &(bdl[x1 * nb + y1]), ldt, &(f77[x1 * lda + y1]), lda);
        }
    }
}

/**
 * @brief Performs Cholesky decomposition using tiled algorithm with OpenMP parallelization
 *
 * This function implements a blocked/tiled version of the DPOTRF (double precision Cholesky
 * factorization) algorithm. It converts the input matrix to a tile layout, performs the
 * factorization in parallel, and converts the result back to the original layout.
 *
 * @details The algorithm uses the following approach:
 * - Allocates temporary memory for tiled matrix storage
 * - Converts input matrix A from column-major to tile layout
 * - Performs parallel Cholesky factorization on tiled data
 * - Converts result back to original column-major layout
 * - Frees temporary memory
 *
 * Memory allocation size calculation:
 * - Block sizes: nb (column blocks), mb (row blocks) are auto-tuned based on problem size
 * - Number of tiles: mt, nt calculated based on matrix dimensions
 * - Total memory includes full tiles plus remainder elements
 *
 * @param A Input/output matrix to be factorized (modified in-place)
 * @param lda Leading dimension of matrix A
 * @param n Pointer to matrix dimension (n x n matrix)
 * @param uplo Specifies upper ('U') or lower ('L') triangular storage
 * @param info Pointer to error information output
 *
 * @note Memory allocation failure results in fallback to reference algorithm(lapack_dpotrf).
 *
 * @note Reference: "A class of parallel tiled linear algebra algorithms for
 *  multicore architectures. Parallel Computing, 35(1), 38-53"
 *  by "Buttari, A., Langou, J., Kurzak, J., & Dongarra, J"
 */
void lapack_dpotrf_var1(char *uplo, aocl_int64_t *n, doublereal *A, aocl_int64_t *lda, aocl_int64_t *info)
{
    *info = 0;
    if(!(*uplo == 'U' || *uplo == 'u' || *uplo == 'L' || *uplo == 'l'))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        aocl_int64_t i__1 = -(*info);
        aocl_blas_xerbla("DPOTRF", &i__1, (ftnlen)6);
        return;
    }
    // Quick return if possible
    if(*n == 0)
    {
        return;
    }

    // Auto-tune parameters based on problem size
    int auto_num_threads;
    aocl_int64_t auto_block_size;

    dpotrf_auto_tune_params(*n, &auto_num_threads, &auto_block_size);

    aocl_int64_t nb = auto_block_size;
    aocl_int64_t mb = auto_block_size;
    aocl_int64_t mt = (*n == 0) ? 0 : (*n - 1) / nb + 1;
    aocl_int64_t nt = (*n == 0) ? 0 : (*n - 1) / nb + 1;
    aocl_int64_t lm1 = *n / mb;
    aocl_int64_t ln1 = *n / nb;
    aocl_int64_t mnt = (ln1 * (1 + lm1)) / 2;
    size_t size = (size_t)(mnt * mb * nb + (*n * (*n % nb))) * sizeof(double);
    double *temp_A = malloc(size);
    A21 = (size_t)(mb * nb) * mnt;
    A12 = (size_t)(mb * nb) * mnt;
    A22 = (size_t)(*n - *n % mb) * (*n % nb) + A12;
    if(temp_A == NULL)
    {
        /*
         * Unable to allocate memory, switching to reference algorithm
         */
        lapack_dpotrf(uplo, n, A, lda, info);
        return;
    }
    char uplo_local = 'L';
    if(*uplo == 'U' || *uplo == 'u')
    {
        uplo_local = 'U';
    }

    if(*n <= PROC_BIND_CLOSE_SIZE || auto_num_threads < PROC_BIND_CLOSE_THREADS)
    {
#pragma omp parallel num_threads(auto_num_threads) proc_bind(close)
#pragma omp single
        {
            // Translate to tile layout.
            matrix_tile(A, *lda, temp_A, nb, mb, &uplo_local, mt, nt, *n, *n, *n);

            // Call to tiled potrf path
            omp_dpotrf(&uplo_local, temp_A, n, lda, mt, mb, nb, *n, *n, info);

            // Get time for matrix_untile
            matrix_untile(A, *lda, temp_A, nb, mb, &uplo_local, mt, nt, *n, *n, *n);
        }
    }
    else
    {
#pragma omp parallel num_threads(auto_num_threads) proc_bind(spread)
#pragma omp single
        {
            // Translate to tile layout.
            matrix_tile(A, *lda, temp_A, nb, mb, &uplo_local, mt, nt, *n, *n, *n);

            // Call to tiled potrf path
            omp_dpotrf(&uplo_local, temp_A, n, lda, mt, mb, nb, *n, *n, info);

            // Get time for matrix_untile
            matrix_untile(A, *lda, temp_A, nb, mb, &uplo_local, mt, nt, *n, *n, *n);
        }
    }
    free(temp_A);
}
#endif
