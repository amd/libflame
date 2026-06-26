/*
    Copyright (C) 2020-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "FLA_f2c.h"
#include "fla_lapack_x86_common.h"
#include "fla_mem.h"

typedef enum dspffrt2_variant_t
{
    DSPFFRT2_UNP_VAR1 = 1,
    DSPFFRT2_UNP_VAR2 = 2,
    DSPFFRT2_UNP_VAR3 = 3,
    DSPFFRT2_UNP_VAR4 = 4,
    DSPFFRT2_UNP_VAR5 = 5,
    DSPFFRT2_UNP_VAR6 = 6
} dspffrt2_variant_t;

typedef struct
{
    dspffrt2_variant_t variant;
    aocl_int64_t onb;
    aocl_int64_t inb;
} dspffrt2_pick_result_t;

static void dspffrt2_fla_def(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                             doublereal *work);
static void dspffrt2_fla_unp_var1(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t nb);
static void dspffrt2_fla_unp_var2(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t nb);
static void dspffrt2_fla_unp_var3(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t nb);
static aocl_int64_t dspffrt2_unp_nb(aocl_int64_t n, aocl_int64_t ncolm);
static void dspffrt2_fla_unp_var4(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t onb, aocl_int64_t inb);
static void dspffrt2_fla_unp_var5(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t onb, aocl_int64_t inb);
static void dspffrt2_fla_unp_var6(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t onb, aocl_int64_t inb);
static dspffrt2_pick_result_t dspffrt2_pick_variant_zen5(aocl_int64_t n, aocl_int64_t ncolm_pc);

extern void DTL_Trace(uint8 ui8LogLevel, uint8 ui8LogType, const int8 *pi8FileName,
                      const int8 *pi8FunctionName, uint32 ui32LineNumber, const int8 *pi8Message);

/*! @brief Partial LDL' factorization without pivoting
    *
    * @details
    * \b Purpose:
    * \verbatim
        DSPFFRT2 computes the partial factorization of a real matrix A
        stored in packed format.
        The factorization has the form
            A = L*D*L**T
        where L is a lower triangular matrix, and D is a diagonal matrix.
        This is an unblocked algorithm.
        The algorthm does not do pivoting and does not handle zero diagonal elements.
        Hence, it may give unexpected outputs for certain inputs.
    \endverbatim

    * @param[in,out] ap
    ap is DOUBLE PRECISION array, dimension (N*(N+1)/2)
    On entry, the lower triangle of the symmetric matrix A, packed columnwise in a
    linear array. The j-th column of A is stored in the array AP as follows:
            AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
    On exit, the block diagonal matrix D and the multipliers used
    to obtain the factor L, stored as a packed triangular matrix overwriting A
    (see below for further details).

    * @param[in] n
    n is integer*. \n
    The order of the matrix A. *n >= 0

    * @param[in] ncolm
    ncolm is integer*. \n
    The number of columns / rows to be factorized. 0 <= *ncolm <= *n

    * @param[in] work
    work is DOUBLE PRECISION array. \n
    Currently an unused buffer

    * @param[in] work2
    work2 is DOUBLE PRECISION array. \n
    Currently an unused buffer

\par Further Details:
   ===================

    * \verbatim
    If input matrix A is of the form,

        ( a  b**T )
    A = (         )
        ( b    C  ), where

    a is the first diagonal element  A, b is a column vector of size n - 1 containing the
    elements from the first column of A excluding the diagonal element,
    C is the lower-right square submatrix of A, and I is the identity matrix,
    then DSPFFRT2 performs ncolm successive factorizations of the form:

        ( a  b**T )   ( a  0 )   ( 1/a       0       )   ( a  b**T )
    A = (         ) = (      ) * (                   ) * (         )
        ( b    C  )   ( b  I )   (  0   C-b*1/a*b**T )   ( 0    I  )

    \endverbatim
    *  */

void dspffrt2_fla(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm, doublereal *work,
                  doublereal *work2)
{
    /* ncolm as fraction of n */
    aocl_int64_t ncolm_pc = (integer)((*ncolm * 100) / *n);

    if(fla_thread_get_num_threads() == 1)
    {
        /* Performance tuned for single-threaded execution */
        if(*n <= FLA_DSPFFRT2_ST__NTHREAD0)
        {
            dspffrt2_fla_def(ap, n, ncolm, work);
        }
        else
        {
            dspffrt2_pick_result_t res = dspffrt2_pick_variant_zen5(*n, ncolm_pc);
            switch(res.variant)
            {
                case DSPFFRT2_UNP_VAR1:
                    dspffrt2_fla_unp_var1(ap, n, ncolm, work, res.inb);
                    break;
                case DSPFFRT2_UNP_VAR2:
                    dspffrt2_fla_unp_var2(ap, n, ncolm, work, res.inb);
                    break;
                case DSPFFRT2_UNP_VAR3:
                    dspffrt2_fla_unp_var3(ap, n, ncolm, work, res.inb);
                    break;
                case DSPFFRT2_UNP_VAR4:
                    dspffrt2_fla_unp_var4(ap, n, ncolm, work, res.onb, res.inb);
                    break;
                case DSPFFRT2_UNP_VAR5:
                    dspffrt2_fla_unp_var5(ap, n, ncolm, work, res.onb, res.inb);
                    break;
                case DSPFFRT2_UNP_VAR6:
                    dspffrt2_fla_unp_var6(ap, n, ncolm, work, res.onb, res.inb);
                    break;
                default:
                    dspffrt2_fla_def(ap, n, ncolm, work);
            }
        }
    }
    else
    {
        if((*n < FLA_SPFFRT2__NTHRESH1)
           || (*n < FLA_SPFFRT2__NTHRESH2 && ncolm_pc < FLA_SPFFRT2__NCOLFRAC_THRESH1)
           || (*ncolm <= FLA_SPFFRT2__NCOLTHRESH))
        {
            /* dspr based implementation for small problem sizes */
            dspffrt2_fla_def(ap, n, ncolm, work);
        }
        else if(ncolm_pc < FLA_SPFFRT2__NCOLFRAC_THRESH2)
        {
            /* Unpacking/packing based variant for smaller ncolm values */
            dspffrt2_fla_unp_var1(ap, n, ncolm, work, 0);
        }
        else
        {
            /* Unpacking/packing based variant for large ncolm values */
            dspffrt2_fla_unp_var2(ap, n, ncolm, work, 0);
        }
    }
}

/*
 * Unpacks selected packed lower-triangular columns into full storage.
 * Writes only the lower part; strict upper entries are left untouched.
 */
void dunpack_fla(doublereal *ap, doublereal *a, aocl_int64_t m, aocl_int64_t start_col,
                 aocl_int64_t end_col, aocl_int64_t lda, aocl_int64_t shift)
{

    aocl_int64_t i;
    aocl_int64_t offset = (start_col * m) - (start_col * (start_col - 1)) / 2;
    doublereal *aptr = ap + offset;

    for(i = start_col; i < end_col; i++)
    {
        aocl_int64_t len = m - i;
        aocl_int64_t ai = i - shift;
        if(len > 0)
        {
            memcpy(&a[ai * lda + ai], aptr, len * sizeof(doublereal));
            aptr += len;
        }
    }
}

/*
 * Packs selected lower-triangular columns from full storage into packed form.
 * Writes only lower data and leaves strict upper entries unused.
 */
void dpack_fla(doublereal *a, doublereal *ap, aocl_int64_t m, aocl_int64_t start_col,
               aocl_int64_t end_col, aocl_int64_t lda, aocl_int64_t shift)
{
    aocl_int64_t i;
    aocl_int64_t offset = (start_col * m) - (start_col * (start_col - 1)) / 2;
    doublereal *aptr = ap + offset;

    for(i = start_col; i < end_col; i++)
    {
        aocl_int64_t len = m - i;
        aocl_int64_t ai = i - shift;
        if(len > 0)
        {
            memcpy(aptr, &a[ai * lda + ai], len * sizeof(doublereal));
            aptr += len;
        }
    }
}

/*
 * LDLT factorization of skinny symmetric matrices (m > n)
 * in unpacked format.
 *
 * Only the lower trapezoidal part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 */
void dsffrk2_fla(doublereal *au, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *lda,
                 doublereal *bt, aocl_int64_t *ldbt)
{
    doublereal d__1;
    aocl_int64_t i__1, i__2, i__3;
    aocl_int64_t k, kc, kcn;
    aocl_int64_t c__1 = 1;
    doublereal r1;

    --au;
    bt -= (1 + *ldbt);
    kc = 1;
    i__3 = *m - *n;
    for(k = 1; k <= *n; k++)
    {
        /* D(k) = -1/A(k,k) */

        /* Skip trailing matrix update if zero diagonal element is encountered */
        if(au[kc] == 0)
            r1 = 0;
        else
            r1 = 1. / au[kc];

        d__1 = -r1;

        i__1 = *n - k;
        i__2 = *m - k;
        kcn = kc + *lda + 1;

        /* Update trailing matrix with rank-1 operation */
        aocl_blas_dger(&i__2, &i__1, &d__1, &au[kc + 1], &c__1, &au[kc + 1], &c__1, &au[kcn], lda);

        /* Compute b**T/a */
#if FLA_ENABLE_AMD_OPT
        fla_dcopy_scal(i__3, d__1, &au[kc + *n - k + 1], 1, &bt[1 + k * *ldbt], 1);
#else
        aocl_blas_dcopy(&i__3, &au[kc + *n - k + 1], &c__1, &bt[1 + k * *ldbt], &c__1);
        aocl_blas_dscal(&i__3, &d__1, &bt[1 + k * *ldbt], &c__1);
#endif

        au[kc] = r1;
        kc = kcn;
    }

    return;
}

/*
 * Baseline packed LDLT path.
 * Uses packed rank-1 updates and serves as the fallback implementation.
 */
void dspffrt2_fla_def(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm, doublereal *work)
{
    doublereal d__1;
    aocl_int64_t i__1, k, kc;
    doublereal r1;

    --ap;

    /* k is the main loop index, increasing from 1 to ncolm in steps of 1 */
    kc = 1;
    for(k = 1; k <= *ncolm; k++)
    {
        /* D(k) = -1/A(k,k) */

        /* Skip trailing matrix update if zero diagonal element is encountered */
        if(ap[kc] == 0)
            r1 = 0;
        else
            r1 = 1. / ap[kc];

        d__1 = -r1;

        /* Update the trailing submatrix */
        /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
        /* A := A - L(k)*D(k)*L(k)**T */
        i__1 = *n - k;
#if FLA_ENABLE_AMD_OPT
        fla_dspr_lower(i__1, d__1, &ap[kc + 1], &ap[kc + *n - k + 1]);
#else
        const aocl_int64_t c__1 = 1;
        aocl_blas_dspr("L", &i__1, &d__1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1]);
#endif

        ap[kc] = r1;
        kc = kc + *n - k + 1;
    }
    return;
}

/*
 * LDLT factorization of symmetric matrices in unpacked format.
 * Blocked algorithm employing GEMMT is used for better
 * performance.
 *
 * Only the lower triangular part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 *
 * Variant 1 does both factorization of (N x ncolm) and
 * trailing matrix update inside the main loop
 */

void dspffrt2_fla_unp_var1(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm, doublereal *work,
                           aocl_int64_t nb)
{
    doublereal d__1 = 1.0;
    aocl_int64_t k, kc;
    aocl_int64_t m;

    doublereal *au, *bt;
    doublereal *mau, *mbt;

    /* Choose block size for the blocked variant */
    nb = nb == 0 ? dspffrt2_unp_nb(*n, *ncolm) : nb;
    nb = (nb > *ncolm) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    aocl_int64_t ldau = FLA_ALIGN(*n, 8);
    aocl_int64_t ldbt = FLA_ALIGN(*n - nb, 8);
    mau = fla_aligned_malloc(ldau * *n * sizeof(doublereal) + ldbt * nb * sizeof(doublereal),
                             FLA_CACHE_LINE_SIZE_BYTES);

    if(mau == NULL)
    {
        /* call default version */
        dspffrt2_fla_def(ap, n, ncolm, work);
        return;
    }

    mbt = mau + ldau * *n;

    dunpack_fla(ap, mau, *n, 0, *n, ldau, 0);

    --ap;
    au = mau - 1;
    bt = mbt - 1;

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 1 to ncolm in steps of nb */
    m = *n;
    for(k = 1; k <= *ncolm; k += nb)
    {
        kc = (k - 1) * ldau + k;
        nb = fla_min(nb, *ncolm - k + 1);

        /* Panel factorization using unblocked variant */
        dsffrk2_fla(&au[kc], &m, &nb, &ldau, &bt[1], &ldbt);
        m -= nb;

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        aocl_blas_dgemm("N", "T", &m, &m, &nb, &d__1, &au[kc + nb], &ldau, &bt[1], &ldbt, &d__1,
                        &au[kc + nb * ldau + nb], &ldau);
#else
        if(m > 0)
            aocl_blas_dgemmt("L", "N", "T", &m, &nb, &d__1, &au[kc + nb], &ldau, &bt[1], &ldbt,
                             &d__1, &au[kc + nb * ldau + nb], &ldau);
#endif

        /* Pack the panel immediately as it will not be updated going on */
        dpack_fla(mau, ap + 1, *n, k - 1, k - 1 + nb, ldau, 0);
    }

    /* Pack the remaining trailing matrix */
    if(*n > *ncolm)
    {
        dpack_fla(mau, ap + 1, *n, *ncolm, *n, ldau, 0);
    }

    fla_aligned_free(mau);
    return;
}

/*
 * LDLT factorization of symmetric matrices in unpacked format.
 * Blocked algorithm employing GEMMT is used for better
 * performance.
 *
 * Only the lower triangular part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 *
 * Variant 2 does factorization of (N x ncolm) in the main loop
 * and the trailing matrix is updated outside the main loop
 */

void dspffrt2_fla_unp_var2(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm, doublereal *work,
                           aocl_int64_t nb)
{
    doublereal d__1 = 1.0;
    aocl_int64_t kc, mg;
    aocl_int64_t k, ni, mp;
    aocl_int64_t kb, ib;

    doublereal *au;
    doublereal *mau;
    doublereal *bt;

    /* Choose block size for the blocked variant */
    nb = nb == 0 ? dspffrt2_unp_nb(*n, *ncolm) : nb;
    nb = (nb > *ncolm) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    aocl_int64_t ldau = FLA_ALIGN(*n, 8);
    aocl_int64_t ldbt = FLA_ALIGN(*n - nb, 8);
    mau = fla_aligned_malloc((ldau + ldbt) * *n * sizeof(doublereal), FLA_CACHE_LINE_SIZE_BYTES);
    if(mau == NULL)
    {
        /* call default version */
        dspffrt2_fla_def(ap, n, ncolm, work);
        return;
    }

    bt = mau + ldau * *n - 1;
    /* Unpack only the panel columns up to ncolm at the start */
    dunpack_fla(ap, mau, *n, 0, *ncolm, ldau, 0);

    --ap;
    au = mau - 1;

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 1 to ncolm in steps of nb */
    mp = *n + nb;
    mg = *n - *ncolm;
    ni = *ncolm;
    for(k = 1; k <= (*ncolm - nb); k += nb)
    {
        mp -= nb;
        kc = (k - 1) * ldau + k;
        ni = ni - nb;
        kb = (k - 1) * ldbt + k;

        /* Panel factorization using unblocked variant */
        dsffrk2_fla(&au[kc], &mp, &nb, &ldau, &bt[kb], &ldbt);

        /* Update trailing matrix within the panel */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        mg = *n - *ncolm + ni;
        aocl_blas_dgemm("N", "T", &mg, &ni, &nb, &d__1, &au[kc + nb], &ldau, &bt[kb], &ldbt, &d__1,
                        &au[kc + nb * ldau + nb], &ldau);
#else
        aocl_blas_dgemmt("L", "N", "T", &ni, &nb, &d__1, &au[kc + nb], &ldau, &bt[kb], &ldbt, &d__1,
                         &au[kc + nb * ldau + nb], &ldau);
        aocl_blas_dgemm("N", "T", &mg, &ni, &nb, &d__1, &au[kc + ni + nb], &ldau, &bt[kb], &ldbt,
                        &d__1, &au[kc + nb * ldau + nb + ni], &ldau);
#endif
        /* Pack the panel immediately as it will not be updated going on */
        dpack_fla(mau, ap + 1, *n, k - 1, k - 1 + nb, ldau, 0);
    }

    /* Process the remaining columns */
    if(k <= *ncolm)
    {
        mp -= nb;
        kc = (k - 1) * ldau + k;
        ib = *ncolm - k + 1;
        kb = (k - 1) * ldbt + k - (nb - ib);

        /* Panel factorization using unblocked variant */
        dsffrk2_fla(&au[kc], &mp, &ib, &ldau, &bt[kb], &ldbt);

        /* Pack the remaining columns of the panel immediately */
        dpack_fla(mau, ap + 1, *n, k - 1, k - 1 + ib, ldau, 0);
    }

    /* Update trailing matrix */
    if(*n > *ncolm)
    {
        aocl_int64_t j, jb, mg_b;
        /* Use a fixed buffer from unused trailing columns space to maximize cache hits */
        doublereal *mau_blk = mau + *ncolm * ldau;

        for(j = *ncolm + 1; j <= *n; j += nb)
        {
            jb = (*n - j + 1 < nb) ? (*n - j + 1) : nb;
            mg_b = *n - j + 1;

            /* Unpack block just when needed into fixed buffer */
            dunpack_fla(ap + 1, mau_blk, *n, j - 1, j - 1 + jb, ldau, j - 1);

#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
            aocl_blas_dgemm("N", "T", &mg_b, &jb, ncolm, &d__1, &au[j], &ldau, &bt[j - nb], &ldbt,
                            &d__1, mau_blk, &ldau);
#else
            aocl_blas_dgemmt("L", "N", "T", &jb, ncolm, &d__1, &au[j], &ldau, &bt[j - nb], &ldbt,
                             &d__1, mau_blk, &ldau);
            if(mg_b > jb)
            {
                aocl_int64_t mg_rem = mg_b - jb;
                aocl_blas_dgemm("N", "T", &mg_rem, &jb, ncolm, &d__1, &au[j + jb], &ldau,
                                &bt[j - nb], &ldbt, &d__1, &mau_blk[jb], &ldau);
            }
#endif

            /* Pack block immediately after its final update from fixed buffer */
            dpack_fla(mau_blk, ap + 1, *n, j - 1, j - 1 + jb, ldau, j - 1);
        }
    }

    fla_aligned_free(mau);
    return;
}

/*
 * Block-size heuristic for unpacked variants 1 and 2.
 * Chooses a tuned block size from matrix size, shape, and threading mode.
 */
static aocl_int64_t dspffrt2_unp_nb(aocl_int64_t n, aocl_int64_t ncolm)
{
    aocl_int64_t nb;

    if(n < FLA_SPFFRT2__BSIZE_NL1)
        nb = FLA_SPFFRT2__BSIZE1;
    else if(n < FLA_SPFFRT2__BSIZE_NL2)
        nb = FLA_SPFFRT2__BSIZE2;
    else
        nb = FLA_SPFFRT2__BSIZE3;

    return nb;
}

/*
 * Blocked unpacked LDLT variant 3 (left-looking).
 * For each block, it first accumulates prior-column contributions,
 * then factorizes the current panel and applies local remainder updates.
 */
void dspffrt2_fla_unp_var3(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm, doublereal *work,
                           aocl_int64_t nb)
{
    doublereal d__1 = 1.0;
    aocl_int64_t k, kc;
    aocl_int64_t m;

    doublereal *bt;
    doublereal *mau;

    /* Choose block size for the blocked variant */
    nb = (nb > *ncolm) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    aocl_int64_t ldau = FLA_ALIGN(*n, 8);
    aocl_int64_t ldbt = FLA_ALIGN(*n - nb, 8);
    mau = fla_aligned_malloc((ldau + ldbt) * *n * sizeof(doublereal), FLA_CACHE_LINE_SIZE_BYTES);

    if(mau == NULL)
    {
        /* call default version */
        dspffrt2_fla_def(ap, n, ncolm, work);
        return;
    }

    bt = mau + ldau * *n;

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 0 to *n in steps of nb */
    for(k = 0; k < *n; k += nb)
    {
        aocl_int64_t ib = fla_min(nb, *n - k);
        kc = k * ldau + k;

        /* Unpack only current nb columns */
        dunpack_fla(ap, mau, *n, k, k + ib, ldau, 0);

        /* Apply all i-1 factors in one go */
        aocl_int64_t d_cols = fla_min(k, *ncolm);
        if(d_cols > 0)
        {
            aocl_int64_t m_rows = *n - k;
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
            aocl_blas_dgemm("N", "T", &m_rows, &ib, &d_cols, &d__1, &mau[k], &ldau, &bt[k - nb],
                            &ldbt, &d__1, &mau[kc], &ldau);
#else
            aocl_blas_dgemmt("L", "N", "T", &ib, &d_cols, &d__1, &mau[k], &ldau, &bt[k - nb], &ldbt,
                             &d__1, &mau[kc], &ldau);
            if(m_rows > ib)
            {
                aocl_int64_t rem_rows = m_rows - ib;
                aocl_blas_dgemm("N", "T", &rem_rows, &ib, &d_cols, &d__1, &mau[k + ib], &ldau,
                                &bt[k - nb], &ldbt, &d__1, &mau[kc + ib], &ldau);
            }
#endif
        }

        /* If within ncolm, factorize the panel */
        if(k < *ncolm)
        {
            aocl_int64_t fib = fla_min(ib, *ncolm - k);
            aocl_int64_t kb = k * ldbt + k - (nb - fib);
            m = *n - k;
            /* Panel factorization using unblocked variant */
            dsffrk2_fla(&mau[kc], &m, &fib, &ldau, &bt[kb], &ldbt);

            /* If fib < ib, update the remaining part of the panel */
            if(fib < ib)
            {
                aocl_int64_t rem_rows = m - fib;
                aocl_int64_t r_cols = ib - fib;
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
                aocl_blas_dgemm("N", "T", &rem_rows, &r_cols, &fib, &d__1, &mau[kc + fib], &ldau,
                                &bt[kb], &ldbt, &d__1, &mau[kc + fib * ldau + fib], &ldau);
#else
                aocl_blas_dgemmt("L", "N", "T", &r_cols, &fib, &d__1, &mau[kc + fib], &ldau,
                                 &bt[kb], &ldbt, &d__1, &mau[kc + fib * ldau + fib], &ldau);
                if(rem_rows > r_cols)
                {
                    aocl_int64_t tr_rows = rem_rows - r_cols;
                    aocl_blas_dgemm("N", "T", &tr_rows, &r_cols, &fib, &d__1, &mau[kc + ib], &ldau,
                                    &bt[kb], &ldbt, &d__1, &mau[kc + fib * ldau + ib], &ldau);
                }
#endif
            }
        }

        /* Pack the panel immediately as it will not be updated going on */
        dpack_fla(mau, ap, *n, k, k + ib, ldau, 0);
    }

    fla_aligned_free(mau);
    return;
}

/*
 * Blocked unpacked LDLT variant 4 (two-level blocking).
 * Uses an outer panel with inner panel factorizations to improve cache reuse,
 * then applies one outer trailing update before packing.
 *
 * The main loop updates the entire trailing matrix after each outer panel factorization.
 */
static void dspffrt2_fla_unp_var4(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t onb, aocl_int64_t inb)
{
    doublereal d__1 = 1.0;
    aocl_int64_t k, kc;
    aocl_int64_t m;

    doublereal *au, *bt;
    doublereal *mau, *mbt;

    aocl_int64_t oib, iib;

    /* Allocate unpacked matrix and do the unpacking */
    aocl_int64_t ldau = FLA_ALIGN(*n, 8);
    aocl_int64_t ldbt = FLA_ALIGN(*n, 8);
    mau = fla_aligned_malloc(ldau * *n * sizeof(doublereal) + ldbt * onb * sizeof(doublereal),
                             FLA_CACHE_LINE_SIZE_BYTES);

    if(mau == NULL)
    {
        /* call default version */
        dspffrt2_fla_def(ap, n, ncolm, work);
        return;
    }

    mbt = mau + ldau * *n;

    dunpack_fla(ap, mau, *n, 0, *n, ldau, 0);

    --ap;
    au = mau - 1;
    bt = mbt - 1;

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 1 to ncolm in steps of nb */
    m = *n;
    for(k = 1; k <= *ncolm; k += onb)
    {
        kc = (k - 1) * ldau + k;
        oib = fla_min(onb, *ncolm - k + 1);

        for(aocl_int64_t jj = 0; jj < oib; jj += inb)
        {
            iib = fla_min(inb, oib - jj);
            aocl_int64_t ikc = (k - 1 + jj) * ldau + (k + jj);
            aocl_int64_t ikbt = jj * ldbt + (jj + iib + 1);
            aocl_int64_t n_right = oib - jj - iib;
            /* Panel factorization using unblocked variant */
            dsffrk2_fla(&au[ikc], &m, &iib, &ldau, &bt[ikbt], &ldbt);
            m -= iib;
            if(n_right > 0)
            {
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
                aocl_blas_dgemm("N", "T", &m, &n_right, &iib, &d__1, &au[ikc + iib], &ldau,
                                &bt[ikbt], &ldbt, &d__1, &au[ikc + iib * ldau + iib], &ldau);
#else
                aocl_blas_dgemmt("L", "N", "T", &n_right, &iib, &d__1, &au[ikc + iib], &ldau,
                                 &bt[ikbt], &ldbt, &d__1, &au[ikc + iib * ldau + iib], &ldau);
                if(m > n_right)
                {
                    aocl_int64_t rem_rows = m - n_right;
                    aocl_blas_dgemm("N", "T", &rem_rows, &n_right, &iib, &d__1,
                                    &au[ikc + iib + n_right], &ldau, &bt[ikbt], &ldbt, &d__1,
                                    &au[ikc + iib * ldau + iib + n_right], &ldau);
                }
#endif
            }
        }

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        aocl_blas_dgemm("N", "T", &m, &m, &oib, &d__1, &au[kc + oib], &ldau, &bt[1 + oib], &ldbt,
                        &d__1, &au[kc + oib * ldau + oib], &ldau);
#else
        aocl_blas_dgemmt("L", "N", "T", &m, &oib, &d__1, &au[kc + oib], &ldau, &bt[1 + oib], &ldbt,
                         &d__1, &au[kc + oib * ldau + oib], &ldau);
#endif
        /* Pack the panel immediately as it will not be updated going on */
        dpack_fla(mau, ap + 1, *n, k - 1, k - 1 + oib, ldau, 0);
    }

    /* Pack the remaining trailing matrix */
    if(*n > *ncolm)
    {
        dpack_fla(mau, ap + 1, *n, *ncolm, *n, ldau, 0);
    }

    fla_aligned_free(mau);
    return;
}

/*
 * Blocked unpacked LDLT variant 5 (two-level blocking).
 * Uses an outer panel with inner panel factorizations to improve cache reuse,
 * then applies one outer trailing update before packing.
 *
 * During factorization, only the panels within the ncolm columns are updated,
 * while the trailing matrix is updated outside the main loop in blocks.
 */
static void dspffrt2_fla_unp_var5(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t onb, aocl_int64_t inb)
{
    doublereal d__1 = 1.0;
    aocl_int64_t k, kc, kb;
    aocl_int64_t m;
    aocl_int64_t mg, ni;

    doublereal *au, *bt;
    doublereal *mau;

    aocl_int64_t oib, iib;

    /* Allocate unpacked matrix and do the unpacking */
    aocl_int64_t ldau = FLA_ALIGN(*n, 8);
    aocl_int64_t ldbt = FLA_ALIGN(*n, 8);
    mau = fla_aligned_malloc(ldau * (*ncolm + onb) * sizeof(doublereal)
                                 + ldbt * *ncolm * sizeof(doublereal),
                             FLA_CACHE_LINE_SIZE_BYTES);

    if(mau == NULL)
    {
        /* call default version */
        dspffrt2_fla_def(ap, n, ncolm, work);
        return;
    }

    bt = mau + ldau * (*ncolm + onb);
    au = mau;

    dunpack_fla(ap, mau, *n, 0, *ncolm, ldau, 0);

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 1 to ncolm in steps of nb */
    m = *n;
    mg = *n - *ncolm;
    ni = *ncolm;
    for(k = 0; k < *ncolm; k += onb)
    {
        kc = k * ldau + k;
        oib = fla_min(onb, *ncolm - k);
        ni -= oib;
        kb = k * ldbt + k;

        for(aocl_int64_t jj = 0; jj < oib; jj += inb)
        {
            iib = fla_min(inb, oib - jj);
            aocl_int64_t ikc = (k + jj) * ldau + (k + jj);
            aocl_int64_t ikbt = (k + jj) * ldbt + (k + jj + iib);
            aocl_int64_t n_right = oib - jj - iib;
            /* Panel factorization using unblocked variant */
            dsffrk2_fla(&au[ikc], &m, &iib, &ldau, &bt[ikbt], &ldbt);
            m -= iib;
            if(n_right > 0)
            {
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
                aocl_blas_dgemm("N", "T", &m, &n_right, &iib, &d__1, &au[ikc + iib], &ldau,
                                &bt[ikbt], &ldbt, &d__1, &au[ikc + iib * ldau + iib], &ldau);
#else
                aocl_blas_dgemmt("L", "N", "T", &n_right, &iib, &d__1, &au[ikc + iib], &ldau,
                                 &bt[ikbt], &ldbt, &d__1, &au[ikc + iib * ldau + iib], &ldau);
                if(m > n_right)
                {
                    aocl_int64_t rem_rows = m - n_right;
                    aocl_blas_dgemm("N", "T", &rem_rows, &n_right, &iib, &d__1,
                                    &au[ikc + iib + n_right], &ldau, &bt[ikbt], &ldbt, &d__1,
                                    &au[ikc + iib * ldau + iib + n_right], &ldau);
                }
#endif
            }
        }

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        mg = *n - *ncolm + ni;
        aocl_blas_dgemm("N", "T", &mg, &ni, &oib, &d__1, &au[kc + oib], &ldau, &bt[kb + oib], &ldbt,
                        &d__1, &au[kc + oib * ldau + oib], &ldau);
#else
        aocl_blas_dgemmt("L", "N", "T", &ni, &oib, &d__1, &au[kc + oib], &ldau, &bt[kb + oib],
                         &ldbt, &d__1, &au[kc + oib * ldau + oib], &ldau);
        aocl_blas_dgemm("N", "T", &mg, &ni, &oib, &d__1, &au[kc + oib + ni], &ldau, &bt[kb + oib],
                        &ldbt, &d__1, &au[kc + oib * ldau + oib + ni], &ldau);
#endif
        /* Pack the panel immediately as it will not be updated going on */
        dpack_fla(mau, ap, *n, k, k + oib, ldau, 0);
    }

    if(*n > *ncolm)
    {
        /* Use a fixed buffer from unused trailing columns space to maximize cache hits */
        doublereal *mau_blk = mau + *ncolm * ldau;

        for(k = *ncolm; k < *n; k += onb)
        {
            oib = fla_min(onb, *n - k);
            mg = *n - k;

            /* Unpack block just when needed into fixed buffer */
            dunpack_fla(ap, mau_blk, *n, k, k + oib, ldau, k);

#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
            aocl_blas_dgemm("N", "T", &mg, &oib, ncolm, &d__1, &au[k], &ldau, &bt[k], &ldbt, &d__1,
                            mau_blk, &ldau);
#else
            aocl_blas_dgemmt("L", "N", "T", &oib, ncolm, &d__1, &au[k], &ldau, &bt[k], &ldbt, &d__1,
                             mau_blk, &ldau);
            if(mg > oib)
            {
                aocl_int64_t mg_rem = mg - oib;
                aocl_blas_dgemm("N", "T", &mg_rem, &oib, ncolm, &d__1, &au[k + oib], &ldau, &bt[k],
                                &ldbt, &d__1, &mau_blk[oib], &ldau);
            }
#endif

            /* Pack block immediately after its final update from fixed buffer */
            dpack_fla(mau_blk, ap, *n, k, k + oib, ldau, k);
        }
    }

    fla_aligned_free(mau);
    return;
}

/*
 * Blocked unpacked LDLT variant 6 (two-level blocking).
 * Uses an outer panel with inner panel factorizations to improve cache reuse,
 * then applies one outer trailing update before packing.
 *
 * During factorization, only the panels within the ncolm columns are updated,
 * while the trailing matrix is updated outside the main loop in one go.
 */
static void dspffrt2_fla_unp_var6(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  doublereal *work, aocl_int64_t onb, aocl_int64_t inb)
{
    doublereal d__1 = 1.0;
    aocl_int64_t k, kc, kb;
    aocl_int64_t m;
    aocl_int64_t mg, ni;

    doublereal *au, *bt;
    doublereal *mau;

    aocl_int64_t oib, iib;

    /* Allocate unpacked matrix and do the unpacking */
    aocl_int64_t ldau = FLA_ALIGN(*n, 8);
    aocl_int64_t ldbt = FLA_ALIGN(*n, 8);
    mau = fla_aligned_malloc(ldau * *n * sizeof(doublereal) + ldbt * *ncolm * sizeof(doublereal),
                             FLA_CACHE_LINE_SIZE_BYTES);

    if(mau == NULL)
    {
        /* call default version */
        dspffrt2_fla_def(ap, n, ncolm, work);
        return;
    }

    bt = mau + ldau * *n;
    au = mau;

    dunpack_fla(ap, mau, *n, 0, *ncolm, ldau, 0);

    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* k is the main loop index, increasing from 1 to ncolm in steps of nb */
    m = *n;
    mg = *n - *ncolm;
    ni = *ncolm;
    for(k = 0; k < *ncolm; k += onb)
    {
        kc = k * ldau + k;
        oib = fla_min(onb, *ncolm - k);
        ni -= oib;
        kb = k * ldbt + k;

        for(aocl_int64_t jj = 0; jj < oib; jj += inb)
        {
            iib = fla_min(inb, oib - jj);
            aocl_int64_t ikc = (k + jj) * ldau + (k + jj);
            aocl_int64_t ikbt = (k + jj) * ldbt + (k + jj + iib);
            aocl_int64_t n_right = oib - jj - iib;
            /* Panel factorization using unblocked variant */
            dsffrk2_fla(&au[ikc], &m, &iib, &ldau, &bt[ikbt], &ldbt);
            m -= iib;
            if(n_right > 0)
            {
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
                aocl_blas_dgemm("N", "T", &m, &n_right, &iib, &d__1, &au[ikc + iib], &ldau,
                                &bt[ikbt], &ldbt, &d__1, &au[ikc + iib * ldau + iib], &ldau);
#else
                aocl_blas_dgemmt("L", "N", "T", &n_right, &iib, &d__1, &au[ikc + iib], &ldau,
                                 &bt[ikbt], &ldbt, &d__1, &au[ikc + iib * ldau + iib], &ldau);
                if(m > n_right)
                {
                    aocl_int64_t rem_rows = m - n_right;
                    aocl_blas_dgemm("N", "T", &rem_rows, &n_right, &iib, &d__1,
                                    &au[ikc + iib + n_right], &ldau, &bt[ikbt], &ldbt, &d__1,
                                    &au[ikc + iib * ldau + iib + n_right], &ldau);
                }
#endif
            }
        }

        /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        mg = *n - *ncolm + ni;
        aocl_blas_dgemm("N", "T", &mg, &ni, &oib, &d__1, &au[kc + oib], &ldau, &bt[kb + oib], &ldbt,
                        &d__1, &au[kc + oib * ldau + oib], &ldau);
#else
        aocl_blas_dgemmt("L", "N", "T", &ni, &oib, &d__1, &au[kc + oib], &ldau, &bt[kb + oib],
                         &ldbt, &d__1, &au[kc + oib * ldau + oib], &ldau);
        aocl_blas_dgemm("N", "T", &mg, &ni, &oib, &d__1, &au[kc + oib + ni], &ldau, &bt[kb + oib],
                        &ldbt, &d__1, &au[kc + oib * ldau + oib + ni], &ldau);
#endif
        /* Pack the panel immediately as it will not be updated going on */
        dpack_fla(mau, ap, *n, k, k + oib, ldau, 0);
    }

    if(*n > *ncolm)
    {
        dunpack_fla(ap, mau, *n, *ncolm, *n, ldau, 0);
        mg = *n - *ncolm;
        kc = *ncolm * ldau + *ncolm;
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        aocl_blas_dgemm("N", "T", &mg, &mg, ncolm, &d__1, &au[*ncolm], &ldau, &bt[*ncolm], &ldbt,
                        &d__1, &au[kc], &ldau);
#else
        aocl_blas_dgemmt("L", "N", "T", &mg, ncolm, &d__1, &au[*ncolm], &ldau, &bt[*ncolm], &ldbt,
                         &d__1, &au[kc], &ldau);
#endif

        dpack_fla(mau, ap, *n, *ncolm, *n, ldau, 0);
    }

    fla_aligned_free(mau);
    return;
}

/*
 * Tuned decision tree for single-threaded execution for Zen5.
 * Chooses the best unpacked variant from matrix size and ncolm percentage.
 */
static dspffrt2_pick_result_t dspffrt2_pick_variant_zen5(aocl_int64_t n, aocl_int64_t ncolm_pc)
{
    if(n <= 1663)
    {
        if(ncolm_pc <= 30)
        {
            if(n <= 659)
            {
                if(n <= 280)
                {
                    if(ncolm_pc <= 18)
                    {
                        if(ncolm_pc <= 0)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                        }
                        else
                        { /* ncolm_pc > 0 */
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 160, 16};
                        }
                    }
                    else
                    { /* ncolm_pc > 18 */
                        if(n <= 276)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                        }
                        else
                        { /* n > 276 */
                            if(ncolm_pc <= 24)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                            else
                            { /* ncolm_pc > 24 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 96, 8};
                            }
                        }
                    }
                }
                else
                { /* n > 280 */
                    if(ncolm_pc <= 18)
                    {
                        if(ncolm_pc <= 17)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                        }
                        else
                        { /* ncolm_pc > 17 */
                            if(n <= 282)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                            }
                            else
                            { /* n > 282 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                        }
                    }
                    else
                    { /* ncolm_pc > 18 */
                        if(ncolm_pc <= 24)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                        }
                        else
                        { /* ncolm_pc > 24 */
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 96, 8};
                        }
                    }
                }
            }
            else
            { /* n > 659 */
                if(ncolm_pc <= 4)
                {
                    if(n <= 968)
                    {
                        if(ncolm_pc <= 3)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                        }
                        else
                        { /* ncolm_pc > 3 */
                            if(n <= 963)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                            else
                            { /* n > 963 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                            }
                        }
                    }
                    else
                    { /* n > 968 */
                        return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 32, 8};
                    }
                }
                else
                { /* ncolm_pc > 4 */
                    if(ncolm_pc <= 8)
                    {
                        if(n <= 966)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                        }
                        else
                        { /* n > 966 */
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                        }
                    }
                    else
                    { /* ncolm_pc > 8 */
                        if(n <= 751)
                        {
                            if(ncolm_pc <= 24)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                            }
                            else
                            { /* ncolm_pc > 24 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 96, 8};
                            }
                        }
                        else
                        { /* n > 751 */
                            if(n <= 1655)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                            }
                            else
                            { /* n > 1655 */
                                if(ncolm_pc <= 18)
                                {
                                    return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                                }
                                else
                                { /* ncolm_pc > 18 */
                                    return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        { /* ncolm_pc > 30 */
            if(n <= 751)
            {
                if(ncolm_pc <= 42)
                {
                    if(n <= 276)
                    {
                        return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                    }
                    else
                    { /* n > 276 */
                        if(n <= 751)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 96, 8};
                        }
                        else
                        { /* n > 751 */
                            if(ncolm_pc <= 31)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                            }
                            else
                            { /* ncolm_pc > 31 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                        }
                    }
                }
                else
                { /* ncolm_pc > 42 */
                    if(n <= 276)
                    {
                        if(ncolm_pc <= 79)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                        }
                        else
                        { /* ncolm_pc > 79 */
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR3, 0, 8};
                        }
                    }
                    else
                    { /* n > 276 */
                        if(n <= 750)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR3, 0, 8};
                        }
                        else
                        { /* n > 750 */
                            if(ncolm_pc <= 80)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR3, 0, 8};
                            }
                            else
                            { /* ncolm_pc > 80 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                        }
                    }
                }
            }
            else
            { /* n > 751 */
                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 96, 8};
            }
        }
    }
    else
    { /* n > 1663 */
        if(ncolm_pc <= 18)
        {
            if(n <= 6095)
            {
                if(n <= 1895)
                {
                    if(ncolm_pc <= 4)
                    {
                        if(n <= 1890)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 32, 8};
                        }
                        else
                        { /* n > 1890 */
                            if(ncolm_pc <= 3)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                            else
                            { /* ncolm_pc > 3 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                            }
                        }
                    }
                    else
                    { /* ncolm_pc > 4 */
                        if(ncolm_pc <= 18)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                        }
                        else
                        { /* ncolm_pc > 18 */
                            if(n <= 1858)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                            }
                            else
                            { /* n > 1858 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                        }
                    }
                }
                else
                { /* n > 1895 */
                    return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                }
            }
            else
            { /* n > 6095 */
                if(ncolm_pc <= 4)
                {
                    if(n <= 13482)
                    {
                        return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                    }
                    else
                    { /* n > 13482 */
                        if(ncolm_pc <= 0)
                        {
                            if(n <= 16396)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR2, 0, 16};
                            }
                            else
                            { /* n > 16396 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                        }
                        else
                        { /* ncolm_pc > 0 */
                            if(n <= 16396)
                            {
                                if(ncolm_pc <= 3)
                                {
                                    return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR6, 96, 8};
                                }
                                else
                                { /* ncolm_pc > 3 */
                                    return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                                }
                            }
                            else
                            { /* n > 16396 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                        }
                    }
                }
                else
                { /* ncolm_pc > 4 */
                    if(n <= 16396)
                    {
                        if(n <= 7989)
                        {
                            if(ncolm_pc <= 7)
                            {
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                            }
                            else
                            { /* ncolm_pc > 7 */
                                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                            }
                        }
                        else
                        { /* n > 7989 */
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                        }
                    }
                    else
                    { /* n > 16396 */
                        return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                    }
                }
            }
        }
        else
        { /* ncolm_pc > 18 */
            if(n <= 16396)
            {
                if(n <= 3525)
                {
                    if(ncolm_pc <= 31)
                    {
                        if(n <= 1858)
                        {
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                        }
                        else
                        { /* n > 1858 */
                            return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
                        }
                    }
                    else
                    { /* ncolm_pc > 31 */
                        return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                    }
                }
                else
                { /* n > 3525 */
                    return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR4, 160, 16};
                }
            }
            else
            { /* n > 16396 */
                return (dspffrt2_pick_result_t){DSPFFRT2_UNP_VAR5, 96, 8};
            }
        }
    }
}

/* dspffrt2_fla */
