/*
    Copyright (c) 2020-2025 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_f2c.h"

extern void z_div(dcomplex *, dcomplex *, dcomplex *);

void zspffrt2_fla_def(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work);
static void zspffrt2_fla_unp_var2(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                                  dcomplex *work);

extern void DTL_Trace(uint8 ui8LogLevel, uint8 ui8LogType, const int8 *pi8FileName,
                      const int8 *pi8FunctionName, uint32 ui32LineNumber, const int8 *pi8Message);

/*! @brief Partial LDL' factorization without pivoting
    *
    * @details
    * \b Purpose:
    * \verbatim
        ZSPFFRT2 computes the partial factorization of a scomplex symmetric matrix A
        stored in packed format.
        The factorization has the form
            A = L*D*L**T
        where L is a lower triangular matrix, and D is a diagonal matrix.
        This is an unblocked algorithm.
        The algorthm does not do pivoting and does not handle zero diagonal elements.
        Hence, it may give unexpected outputs for certain inputs.
    \endverbatim

    * @param[in,out] ap
    ap is COMPLEX*16 array, dimension (N*(N+1)/2)
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
    work is COMPLEX*16 array. \n
    Currently an unused buffer

    * @param[in] work2
    work2 is COMPLEX*16 array. \n
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
    then ZSPFFRT2 performs ncolm successive factorizations of the form:

        ( a  b**T )   ( a  0 )   ( 1/a       0       )   ( a  b**T )
    A = (         ) = (      ) * (                   ) * (         )
        ( b    C  )   ( b  I )   (  0   C-b*1/a*b**T )   ( 0    I  )

    \endverbatim
    *  */
void zspffrt2_fla(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work,
                  dcomplex *work2)
{
    /* ncolm as fraction of n */
    aocl_int64_t ncolm_pc = (integer)((*ncolm * 100) / *n);
    if((*n > (FLA_SPFFRT2__NTHRESH1 - 1)) && (ncolm_pc >= FLA_SPFFRT2__NCOLFRAC_THRESH3))
    {
        /* Unpacking/packing based variant for small n & ncolm values */
        zspffrt2_fla_unp_var2(ap, n, ncolm, work);
    }
    else if(*n > FLA_SPFFRT2__NTHRESH3)
    {
        /* Unpacking/packing based variant for large n values */
        zspffrt2_fla_unp_var2(ap, n, ncolm, work);
    }
    else
    {
        /* zspr based implementation for small problem sizes */
        zspffrt2_fla_def(ap, n, ncolm, work);
    }
    return;
}

/*
 *  Unpacking the input packed symmetric matrix into lower
 *  triangular part of unpacked full matrix.
 *  The strictly upper triangular part is left untouched.
 */
void zunpack_fla(dcomplex *ap, dcomplex *a, aocl_int64_t m, aocl_int64_t n,
                 aocl_int64_t lda)
{
    aocl_int64_t i, j;
    dcomplex *aptr = ap;

    for(i = 0; i < n; i++)
    {
        for(j = i; j < m; j++)
        {
            a[i * lda + j] = *aptr++;
        }
    }

    return;
}

/*
 *  Packing the lower triangular part of input symmetric matrix into
 *  packed matrix.
 *  The strictly upper triangular parts of the input and output are
 *  left unused and untouched respectiely.
 */
void zpack_fla(dcomplex *a, dcomplex *ap, aocl_int64_t m, aocl_int64_t n,
               aocl_int64_t lda)
{
    aocl_int64_t i, j;
    dcomplex *aptr = ap;

    for(i = 0; i < n; i++)
    {
        for(j = i; j < m; j++)
        {
            *aptr++ = a[i * lda + j];
        }
    }

    return;
}

/*
 * LDLT factorization of skinny symmetric matrices (M > N)
 * in unpacked format.
 *
 * Only the lower trapezoidal part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 */
void zsffrk2_fla(dcomplex *au, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *lda,
                 dcomplex *bt, aocl_int64_t *ldbt)
{
    dcomplex z__1;
    aocl_int64_t i__1, i__2, i__3;
    aocl_int64_t k, kc, kcn;
    aocl_int64_t c__1 = 1;
    dcomplex r1;
    dcomplex c_b1 = {1., 0.};

    --au;
    --bt;
    kc = 1;
    i__3 = *m - *n;
    for(k = 1; k <= *n; k++)
    {
        /* D(k) = -1/A(k,k) */

        /* Skip trailing matrix update if zero diagonal element is encountered */
        if(au[kc].real == 0 && au[kc].imag == 0)
        {
            z__1.real = 0;
            z__1.imag = 0;
        }
        else
        {
            z_div(&z__1, &c_b1, &au[kc]);
        }

        r1.real = -z__1.real;
        r1.imag = -z__1.imag;

        i__1 = *n - k;
        i__2 = *m - k;
        kcn = kc + *lda + 1;

        /* Update trailing matrix with rank-1 operation */
        aocl_blas_zgeru(&i__2, &i__1, &r1, &au[kc + 1], &c__1, &au[kc + 1], &c__1, &au[kcn], lda);

        /* Compute b**T/a */
        aocl_blas_zcopy(&i__3, &au[kc + *n - k + 1], &c__1, &bt[k], ldbt);
        aocl_blas_zscal(&i__3, &r1, &bt[k], ldbt);

        au[kc].real = z__1.real;
        au[kc].imag = z__1.imag;
        kc = kcn;
    }

    return;
}

/*
 * LDLT factorization of symmetric matrices (M > N)
 * in unpacked format.
 * Blocked algorithm employing GEMMT is used for better
 * performance.
 *
 * Only the lower triangular part of the matrix is updated.
 * The strictly upper triangular part is left untouched.
 *
 * Variant 2 does factorization of (N x ncolm) in the main loop
 * and the trailing matrix is updated outside the main loop
 */
void zspffrt2_fla_unp_var2(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm,
                           dcomplex *work)
{
    dcomplex d__1 = {1., 0.};
    aocl_int64_t kc, mg, nb;
    aocl_int64_t k, ni, mp;

    dcomplex *au;
    dcomplex *mau;

    /* Choose block size for the blocked variant */
    if(*n < FLA_SPFFRT2__BSIZE_NL1)
        nb = FLA_SPFFRT2__BSIZE1;
    else if(*n < FLA_SPFFRT2__BSIZE_NL2)
        nb = FLA_SPFFRT2__BSIZE2;
    else
        nb = FLA_SPFFRT2__BSIZE3;
    nb = (nb > *ncolm) ? *ncolm : nb;

    /* Allocate unpacked matrix and do the unpacking */
    mau = NULL;
    mau = (dcomplex *) malloc(*n * *n * sizeof(dcomplex));
    if(mau == NULL)
    {
        /* call default version */
        zspffrt2_fla_def(ap, n, ncolm, work);
        return;
    }

    zunpack_fla(ap, mau, *n, *n, *n);

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
        kc = k * *n - *n + k;
        ni = ni - nb;

        /* Panel factorization using unblocked variant */
        zsffrk2_fla(&au[kc], &mp, &nb, n, &au[kc + nb * *n], n);

        /* Update trailing matrix within the panel */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
        mg = *n - *ncolm + ni;
        aocl_blas_zgemm("N", "N", &mg, &ni, &nb, &d__1, &au[kc + nb], n, &au[kc + nb * *n], n,
                        &d__1,
               &au[kc + nb * *n + nb], n);
#else
        aocl_blas_zgemmt("L", "N", "N", &ni, &nb, &d__1, &au[kc + nb], n, &au[kc + nb * *n], n, &d__1,
                &au[kc + nb * *n + nb], n);
        aocl_blas_zgemm("N", "N", &mg, &ni, &nb, &d__1, &au[kc + ni + nb], n, &au[kc + nb * *n], n,
                &d__1, &au[kc + nb * *n + nb + ni], n);
#endif
    }

    /* Process the remaining columns */
    if(k <= *ncolm)
    {
        mp -= nb;
        kc = k * *n - *n + k;
        nb = *ncolm - k + 1;

        /* Panel factorization using unblocked variant */
        zsffrk2_fla(&au[kc], &mp, &nb, n, &au[kc + nb * *n], n);
    }

    /* Update trailing matrix */
#ifndef FLA_ENABLE_BLAS_EXT_GEMMT
    mg = *n - *ncolm;
    aocl_blas_zgemm("N", "N", &mg, &mg, ncolm, &d__1, &au[*ncolm + 1], n, &au[*ncolm * *n + 1], n, &d__1,
           &au[*ncolm + *ncolm * *n + 1], n);
#else
    aocl_blas_zgemmt("L", "N", "N", &mg, ncolm, &d__1, &au[*ncolm + 1], n, &au[*ncolm * *n + 1], n, &d__1,
            &au[*ncolm + *ncolm * *n + 1], n);
#endif

    /* Pack the updated matrix */
    zpack_fla(mau, &ap[1], *n, *n, *n);

    free(mau);
    return;
}

/*
 * AXPY based optimized version of rank-1 update for packedm
 * matrices.
 */

int lzspr_(char *uplo, aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx,
           dcomplex *ap, dcomplex *work)
{
    aocl_int64_t incw = 1;
    aocl_int64_t k, kn, nz;

    ap--;
    work--;
    x--;

    for(k = 1; k <= *n; k++)
    {
        work[k].real = alpha->real * x[k].real - alpha->imag * x[k].imag;
        work[k].imag = alpha->real * x[k].imag + alpha->imag * x[k].real;
    }

    kn = 1;
    for(k = 1; k <= *n; k++)
    {
        nz = *n - k + 1;
        aocl_blas_zaxpy(&nz, &work[k], &x[k], &incw, &ap[kn], &incw);
        kn = kn + *n - k + 1;
    }

    return 0;
}

void zspffrt2_fla_def(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work)
{
    dcomplex z__1;
    aocl_int64_t i__1, k, kc;
    aocl_int64_t c__1 = 1;
    dcomplex r1;
    dcomplex c_b1 = {1., 0.};

    --ap;
    /* Factorize A as L*D*L**T using the lower triangle of A */
    /* K is the main loop index, increasing from 1 to ncolm in steps of 1 */
    kc = 1;
    for(k = 1; k <= *ncolm; k++)
    {
        /* Update the trailing submatrix */
        /* W(k) = L(k)*D(k) */
        /* where L(k) is the k-th column of L */

        /* Skip trailing matrix update if zero diagonal element is encountered */
        if(ap[kc].real == 0 && ap[kc].imag == 0)
        {
            z__1.real = 0;
            z__1.imag = 0;
        }
        else
        {
            z_div(&z__1, &c_b1, &ap[kc]);
        }

        r1.real = -z__1.real;
        r1.imag = -z__1.imag;

        /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
        /* A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */
        i__1 = *n - k;
        lzspr_("Lower", &i__1, &r1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1], work);

        ap[kc].real = z__1.real;
        ap[kc].imag = z__1.imag;

        kc = kc + *n - k + 1;
    }
    return;
}
/* zspffrt2_fla */
