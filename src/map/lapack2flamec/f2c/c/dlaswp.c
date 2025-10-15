/* ../netlib/dlaswp.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/*
*     Modifications Copyright (c) 2025 Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "FLA_f2c.h" /* > \brief \b DLASWP performs a series of row interchanges on a general rectangular matrix. */

#define DLASWP_LARGE_THRESHOLD_MIN 1024
#define DLASWP_TILE_SIZE_SMALL 32
#define DLASWP_TILE_SIZE_LARGE 64
#define DLASWP_TILE_MASK_SMALL 31  /* DLASWP_TILE_SIZE_SMALL - 1 */
#define DLASWP_TILE_MASK_LARGE 63  /* DLASWP_TILE_SIZE_LARGE - 1 */

/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLASWP + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, K1, K2, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASWP performs a series of row interchanges on the matrix A. */
/* > One row interchange is initiated for each of rows K1 through K2 of A. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the matrix of column dimension N to which the row */
/* > interchanges will be applied. */
/* > On exit, the permuted matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > \endverbatim */
/* > */
/* > \param[in] K1 */
/* > \verbatim */
/* > K1 is INTEGER */
/* > The first element of IPIV for which a row interchange will */
/* > be done. */
/* > \endverbatim */
/* > */
/* > \param[in] K2 */
/* > \verbatim */
/* > K2 is INTEGER */
/* > (K2-K1+1) is the number of elements of IPIV for which a row */
/* > interchange will be done. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (K1+(K2-K1)*f2c_abs(INCX)) */
/* > The vector of pivot indices. Only the elements in positions */
/* > K1 through K1+(K2-K1)*f2c_abs(INCX) of IPIV are accessed. */
/* > IPIV(K1+(K-K1)*f2c_abs(INCX)) = L implies rows K and L are to be */
/* > interchanged. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of IPIV. If INCX */
/* > is negative, the pivots are applied in reverse order. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup doubleOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Modified by */
/* > R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dlaswp_(integer *n, doublereal *a, integer *lda, integer *k1, integer *k2, integer *ipiv,
             integer *incx)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlaswp inputs: n %" FLA_IS ", lda %" FLA_IS ", k1 %" FLA_IS ", k2 %" FLA_IS
                      ", incx %" FLA_IS "",
                      *n, *lda, *k1, *k2, *incx);
    /* System generated locals */
    integer a_dim1, a_offset, i__4;
    /* Local variables */
    integer i__, j, k, i1, i2, n32, n64, ip, ix, ix0, inc;
    doublereal *colptr, temp;
    /* -- LAPACK auxiliary routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Interchange row I with row IPIV(K1+(I-K1)*f2c_abs(INCX)) for each of rows */
    /* K1 through K2. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    /* Function Body */
    if(*incx > 0)
    {
        ix0 = *k1;
        i1 = *k1;
        i2 = *k2;
        inc = 1;
    }
    else if(*incx < 0)
    {
        ix0 = *k1 + (*k1 - *k2) * *incx;
        i1 = *k2;
        i2 = *k1;
        inc = -1;
    }
    else
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }

    /* Process in tiles if n < 1024 */
    if(*n < DLASWP_LARGE_THRESHOLD_MIN)
    {
        n32 = *n & ~DLASWP_TILE_MASK_SMALL;  /* Equivalent to (*n / 32) * 32 but faster */
        if(n32 != 0)
        {
            for(j = 1; j <= n32; j += DLASWP_TILE_SIZE_SMALL)
            {
                ix = ix0;
                for(i__ = i1; inc < 0 ? i__ >= i2 : i__ <= i2; i__ += inc)
                {
                    ip = ipiv[ix];
                    if(ip != i__)
                    {
                        i__4 = j + DLASWP_TILE_MASK_SMALL;
                        for(k = j; k <= i__4; ++k)
                        {
                            colptr = &a[k * a_dim1];
                            temp = colptr[i__];
                            colptr[i__] = colptr[ip];
                            colptr[ip] = temp;
                            /* L10: */
                        }
                    }
                    ix += *incx;
                    /* L20: */
                }
                /* L30: */
            }
        }
        if(n32 != *n)
        {
            ++n32;
            ix = ix0;
            for(i__ = i1; inc < 0 ? i__ >= i2 : i__ <= i2; i__ += inc)
            {
                ip = ipiv[ix];
                if(ip != i__)
                {
                    for(k = n32; k <= *n; ++k)
                    {
                        colptr = &a[k * a_dim1];
                        temp = colptr[i__];
                        colptr[i__] = colptr[ip];
                        colptr[ip] = temp;
                        /* L40: */
                    }
                }
                ix += *incx;
                /* L50: */
            }
        }
    }
    else
    {
        /* Process columns in tiles to improve cache locality. Use a larger tile
           size on modern CPUs. */
        n64 = *n & ~DLASWP_TILE_MASK_LARGE;  /* Equivalent to (*n / 64) * 64 but faster */
#ifdef FLA_OPENMP_MULTITHREADING
/*
 #pragma omp parallel for:
    This is a combined directive that creates a team of threads and then divides the iterations of
 the following for loop among them. schedule(static): With static scheduling, the iterations are
 divided into approximately equal-sized contiguous blocks, and each block is assigned to a thread.
 if (condition):
    If the condition true, it executed in parallel. If the condition is zero (false), the for loop
 will be executed sequentially by a single thread. */
#pragma omp parallel for schedule(static) private(ix, i__, ip, k, i__4, colptr, temp) firstprivate(ix0, i1, i2, inc, a_dim1)
#endif
        for(j = 1; j <= n64; j += DLASWP_TILE_SIZE_LARGE)
        {
            ix = ix0;
            for(i__ = i1; inc < 0 ? i__ >= i2 : i__ <= i2; i__ += inc)
            {
                ip = ipiv[ix];
                if(ip != i__)
                {
                    i__4 = j + DLASWP_TILE_MASK_LARGE;
                    for(k = j; k <= i__4; ++k)
                    {
                        colptr = &a[k * a_dim1];
                        temp = colptr[i__];
                        colptr[i__] = colptr[ip];
                        colptr[ip] = temp;
                        /* L10: */
                    }
                }
                ix += *incx;
                /* L20: */
            }
            /* L30: */
        }
        if (n64 < *n)
        {
            ix = ix0;
            n64++;
            for(i__ = i1; inc < 0 ? i__ >= i2 : i__ <= i2; i__ += inc)
            {
                ip = ipiv[ix];
                if(ip != i__)
                {
                    for(k = n64; k <= *n; ++k)
                    {
                        colptr = &a[k * a_dim1];
                        temp = colptr[i__];
                        colptr[i__] = colptr[ip];
                        colptr[ip] = temp;
                        /* L40: */
                    }
                }
                ix += *incx;
                /* L50: */
            }
        }
    }

    AOCL_DTL_TRACE_LOG_EXIT
    /* End of DLASWP */
}
/* dlaswp_ */
