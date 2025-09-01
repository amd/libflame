/* dgeqrf.f -- translated by f2c (version 20000121). You must link the resulting object file with
 * the libraries: -lf2c -lm (in that order) */
/******************************************************************************
 * Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#include "FLA_f2c.h" /* Table of constant values */
#include "fla_lapack_x86_common.h"
#if !FLA_ENABLE_AMD_OPT
static integer c__1 = 1;
#endif
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

extern int fla_thread_get_num_threads();

integer get_block_size_dgeqrf(integer *m, integer *n);

#if FLA_ENABLE_AMD_OPT
static void dgeqrf_mt_large(integer gm, integer gn, doublereal *a, integer lda, doublereal *tau,
                            doublereal *work, integer nthreads, integer *info);
static integer dgeqrf_mt_large_num_threads(integer gm, integer gn);
static integer dgeqrf_mt_large_lwork(integer gm, integer gn, integer num_threads);
void dgeqr2_fla(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                integer *);
void dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *,
             integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
void dlarft_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
             doublereal *, integer *);
#endif

/* > \brief \b DGEQRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEQRF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQRF computes a QR factorization of a real M-by-N matrix A: */
/* > */
/* > A = Q * ( R ), */
/* > ( 0 ) */
/* > */
/* > where: */
/* > */
/* > Q is a M-by-M orthogonal matrix;
 */
/* > R is an upper-triangular N-by-N matrix;
 */
/* > 0 is a (M-N)-by-N zero matrix, if M > N. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the fla_min(M,N)-by-N upper trapezoidal matrix R (R is */
/* > upper triangular if m >= n);
the elements below the diagonal, */
/* > with the array TAU, represent the orthogonal matrix Q as a */
/* > product of fla_min(m,n) elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
/* > For optimum performance LWORK >= N*NB, where NB is */
/* > the optimal blocksize. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup doubleGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k), where k = fla_min(m,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar, and v is a real vector with */
/* > v(1:i-1) = 0 and v(i) = 1;
v(i+1:m) is stored on exit in A(i+1:m,i), */
/* > and tau in TAU(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */

void dgeqrf_fla(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau,
                doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, k, nbmin, iinfo;
    extern /* Subroutine */
        void
        dgeqr2_fla(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                   integer *);
    integer ib, nb;
    extern /* Subroutine */
        void
        dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *,
                integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *);
    integer nx;
    extern /* Subroutine */
        void
        dlarft_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ldwork, lwkopt;
    logical lquery;
    integer iws;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

#if AOCL_FLA_PROGRESS_H
    AOCL_FLA_PROGRESS_VAR;
#endif
    --tau;
    --work;
    /* Function Body */
    *info = 0;
#if FLA_ENABLE_AMD_OPT
    nb = get_block_size_dgeqrf(m, n);
#else
    nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1);
#endif

    lwkopt = *n * nb;
#if FLA_ENABLE_AMD_OPT && FLA_OPENMP_MULTITHREADING
    /* If the matrix is large, use the multithreaded version */
    /* Need to tune this threshold */
    integer opt_mt_threads = 1;
    if(*m >= FLA_DGEQRF_MT_LARGE_M_THRESH && *n >= FLA_DGEQRF_MT_LARGE_N_THRESH)
    {
        opt_mt_threads = dgeqrf_mt_large_num_threads(*m, *n);
        lwkopt = dgeqrf_mt_large_lwork(*m, *n, opt_mt_threads);
    }
#endif
    work[1] = (doublereal)lwkopt;
    lquery = *lwork == -1;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    else if(*lwork < fla_max(1, *n) && !lquery)
    {
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGEQRF", &i__1, (ftnlen)6);
        return;
    }
    else if(lquery)
    {
        return;
    }
    /* Quick return if possible */
    k = fla_min(*m, *n);
    if(k == 0)
    {
        work[1] = 1.;
        return;
    }
    nbmin = 2;
    nx = 0;
    iws = *n;
#if AOCL_FLA_PROGRESS_H
    progress_step_count = 0;
#ifndef FLA_ENABLE_WINDOWS_BUILD
    if(!aocl_fla_progress_ptr)
        aocl_fla_progress_ptr = aocl_fla_progress;
#endif
#endif

/* Path for small sizes */
#if FLA_ENABLE_AMD_OPT
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2) && *m <= FLA_GEQRF_STHRESH && *n <= FLA_GEQRF_STHRESH)
    {
        fla_dgeqrf_small(m, n, &a[a_offset], lda, &tau[1], &work[1]);
        return;
    }
#if FLA_OPENMP_MULTITHREADING
    else if(*m >= FLA_DGEQRF_MT_LARGE_M_THRESH && *n >= FLA_DGEQRF_MT_LARGE_N_THRESH
            && *lwork >= lwkopt)
    {
        dgeqrf_mt_large(*m, *n, &a[a_offset], *lda, &tau[1], &work[1], opt_mt_threads, info);
        return;
    }
#endif
#endif

    if(nb > 1 && nb < k)
    {
        /* Determine when to cross over from blocked to unblocked code. */
        /* Computing MAX */
        i__1 = 0;
        i__2 = ilaenv_(&c__3, "DGEQRF", " ", m, n, &c_n1, &c_n1); // , expr subst
        nx = fla_max(i__1, i__2);
        if(nx < k)
        {
            /* Determine if workspace is large enough for blocked code. */
            ldwork = *n;
            iws = ldwork * nb;
            if(*lwork < iws)
            {
                /* Not enough workspace to use optimal NB: reduce NB and */
                /* determine the minimum value of NB. */
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = ilaenv_(&c__2, "DGEQRF", " ", m, n, &c_n1, &c_n1); // , expr subst
                nbmin = fla_max(i__1, i__2);
            }
        }
    }
    if(nb >= nbmin && nb < k && nx < k)
    {
        /* Use blocked code initially */
        i__1 = k - nx;
        i__2 = nb;
        for(i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            /* Computing MIN */
            i__3 = k - i__ + 1;
            ib = fla_min(i__3, nb);
            /* Compute the QR factorization of the current block */
            /* A(i:m,i:i+ib-1) */
            i__3 = *m - i__ + 1;
#if AOCL_FLA_PROGRESS_H

            if(aocl_fla_progress_ptr)
            {
                progress_step_count += ib;
                AOCL_FLA_PROGRESS_FUNC_PTR("DGEQRF", 6, &progress_step_count, &progress_thread_id,
                                           &progress_total_threads);
            }

#endif

            dgeqr2_fla(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
            if(i__ + ib <= *n)
            {
                /* Form the triangular factor of the block reflector */
                /* H = H(i) H(i+1) . . . H(i+ib-1) */

                i__3 = *m - i__ + 1;
                dlarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__],
                        &work[1], &ldwork);
                /* Apply H**T to A(i:m,i+ib:n) from the left */
                i__3 = *m - i__ + 1;
                i__4 = *n - i__ - ib + 1;
                dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &i__4, &ib,
                        &a[i__ + i__ * a_dim1], lda, &work[1], &ldwork,
                        &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib + 1], &ldwork);
            }
            /* L10: */
        }
    }
    else
    {
        i__ = 1;
    }
    /* Use unblocked code to factor the last or only block. */
    if(i__ <= k)
    {
        i__2 = *m - i__ + 1;
        i__1 = *n - i__ + 1;
#if AOCL_FLA_PROGRESS_H

        if(aocl_fla_progress_ptr)
        {
            progress_step_count = k;
            AOCL_FLA_PROGRESS_FUNC_PTR("DGEQRF", 6, &progress_step_count, &progress_thread_id,
                                       &progress_total_threads);
        }

#endif

        dgeqr2_fla(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
    }
    work[1] = (doublereal)iws;
    return;
    /* End of DGEQRF */
}

integer get_block_size_dgeqrf(integer *m, integer *n)
{
    integer block_size;

    /* Set block_size=32 for small sizes */
    if(*m <= 17 && *n <= 17)
    {
        block_size = 32;
    }
    else
    {
        int num_threads = fla_thread_get_num_threads();

        if(num_threads == 1)
        {
            if(*m >= 2000 && *n >= 2000)
            {
                block_size = 64;
            }
            else if(*m >= 1000 && *n >= 1000)
            {
                block_size = 60;
            }
            else
            {
                block_size = 32;
            }
        }
        else
        {
            block_size = 32;
        }
    }

    return block_size;
}

/* dgeqrf_ */

#if FLA_ENABLE_AMD_OPT && FLA_OPENMP_MULTITHREADING

/* Returns the optimial number of threads to use for dgeqrf_mt_large */
static integer dgeqrf_mt_large_num_threads(integer gm, integer gn)
{
    long long num_elems = (long long)gm * gn;
    integer opt_nthreads;
    if(num_elems <= FLA_DGEQRF_MT_THRESHOLD_8_THREADS)
    {
        opt_nthreads = 8;
    }
    else
    {
        opt_nthreads = 64;
    }
    return fla_min(fla_thread_get_num_threads(), opt_nthreads);
}

/* Returns the number of containers needed for given nb items of size item_size
 * to fit in a container of size container_size */
static inline size_t containers_needed_for_size(integer nb, size_t item_size, size_t container_size)
{
    size_t size_needed = item_size * nb;
    return (size_needed - 1) / container_size + 1;
}

/* Returns the optimal lwork for multithreaded dgeqrf */
static integer dgeqrf_mt_large_lwork(integer gm, integer gn, integer num_threads)
{

    integer nb = FLA_DGEQRF_MT_LARGE_PANEL_SIZE;

    /* Number of panels */
    integer nt = ((gn - 1) / nb) + 1;

    /* Storage to keep triangular factors for each panel */
    integer triangular_factor_req = (nb * nb) * nt;

    integer per_thread_work_req = (nb * nb) * num_threads;

    /* Dependency list to keep track of which panels are ready */
    integer depend_list_req = containers_needed_for_size(nt, sizeof(uint8_t), sizeof(doublereal));

    integer lwork = triangular_factor_req + per_thread_work_req + depend_list_req;

    return lwork;
}

static void dgeqrf_large_mt_thread_fn(integer m, integer n, doublereal *a, integer lda,
                                      doublereal *tau, doublereal *T_storage, uint8_t *a_dependency,
                                      integer nb, integer k, doublereal *per_thread_work,
                                      integer per_thread_work_size, integer *info)
{
    /* Get the current thread number */
    integer thread_num = omp_get_thread_num();
    integer iinfo;
    /* Get the thread-local workspace */
    doublereal *t_work = per_thread_work + thread_num * per_thread_work_size;

    integer ldwork = nb;

    /* Apply block reflector to all panels in left
     * columns of the current block.
     */
    for(integer i = 0; i < fla_min(m, k); i += nb)
    {
        /* If the block reflector is not ready, wait for it */
        while(a_dependency[i / nb] == 0)
        {
            /* Wait for all memory operation to be completed */
            __asm__ __volatile__("lfence");
            /* Use pause instruction for busy waiting */
            __asm__ __volatile__("pause");
        }

        /* Apply the block reflector to the current panel */
        integer mvai = m - i;
        integer n_reflector = fla_min(mvai, nb);
        dlarfb_("Left", "Transpose", "Forward", "Columnwise", &mvai, &n, &n_reflector,
                &a[i + i * lda], &lda, T_storage + (nb * i), &nb, &a[k * lda + i], &lda, t_work,
                &ldwork);
    }

    integer mvak = m - k;

    /* If there are no more rows to process, return */
    if(mvak <= 0)
    {
        return;
    }

    /* Factor the current panel */
    dgeqr2_fla(&mvak, &n, &a[k * lda + k], &lda, &tau[k], t_work, &iinfo);

    /* If the factorization failed, set the info and return */
    if(iinfo != 0)
    {
        *info = iinfo;
        return;
    }

    /* Generate the triangular factor of the block reflector */
    integer k_reflector = fla_min(mvak, n);
    dlarft_("Forward", "Columnwise", &mvak, &k_reflector, &a[k * lda + k], &lda, &tau[k],
            T_storage + (nb * k), &nb);

    /* Set the block reflector as ready */
    a_dependency[k / nb] = 1;
    /* Ensure all memory operations are completed */
    __asm__ __volatile__("sfence");
}

static void dgeqrf_mt_large(integer gm, integer gn, doublereal *a, integer lda, doublereal *tau,
                            doublereal *work, integer nthreads, integer *info)
{

    integer nb = FLA_DGEQRF_MT_LARGE_PANEL_SIZE;

    integer iinfo = 0;

    /* Number of panels */
    integer nt = ((gn - 1) / nb) + 1;

    /* Storage to keep triangular factors for each panel */
    doublereal *T_storage = work;

    /* Dependency list to keep track of which panels are ready */
    uint8_t *a_dependency = (uint8_t *)(T_storage + (nb * nb) * nt);

    /* Per-thread workspace */
    doublereal *per_thread_work
        = ((doublereal *)a_dependency)
          + containers_needed_for_size(nt, sizeof(uint8_t), sizeof(doublereal));

    /* Initialize the dependency list to 0 */
    memset(a_dependency, 0, nt);

    /* Work size for each thread */
    integer per_thread_work_size = (nb * nb);

/* Schedule threads to process panels
 * Panels are assigned to threads in a round-robin manner
 */
#pragma omp parallel for num_threads(nthreads) schedule(static, 1) proc_bind(close)
    for(integer i = 0; i < gn; i += nb)
    {
        /* Get the number of coulmns in current panel */
        integer n_thread = fla_min(nb, gn - i);
        dgeqrf_large_mt_thread_fn(gm, n_thread, a, lda, tau, T_storage, a_dependency, nb, i,
                                  per_thread_work, per_thread_work_size, &iinfo);
    }

    /* Set the info */
    *info = iinfo;
}

#endif
