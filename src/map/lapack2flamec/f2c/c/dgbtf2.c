/* ../netlib/dgbtf2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
   Modifications Copyright (c) 2024 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h" /* Table of constant values */
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

static integer c__1 = 1;
static doublereal c_b9 = -1.;

#if FLA_ENABLE_AOCL_BLAS

/* This function is an implementation of dgbtf2 using AOCL-BLAS compute kernels */
void dgbtf2_aocl_blas_ver(integer *m, integer *n, integer *kl, integer *ku, doublereal *ab,
                          integer *ldab, integer *ipiv, integer *info)
{
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    integer i__, j, km, jp, ju, kv;
    doublereal alpha;
    doublereal *x, *y, *r;
    dim_t index;

    /* Make a copy of AOCL-BLAS framework context. This information is needed to query the
     * architecture specific details of compute kernel */
    cntx_t *cntx = bli_gks_query_cntx();

    /* Query names of compute kernel from AOCL-BLAS framework context */
    damaxv_ker_ft idamax_blas_ptr = bli_cntx_get_l1v_ker_dt(BLIS_DOUBLE, BLIS_AMAXV_KER, cntx);
    dswapv_ker_ft dswap_blas_ptr = bli_cntx_get_l1v_ker_dt(BLIS_DOUBLE, BLIS_SWAPV_KER, cntx);
    dscalv_ker_ft dscal_blas_ptr = bli_cntx_get_l1v_ker_dt(BLIS_DOUBLE, BLIS_SCALV_KER, cntx);
    daxpyv_ker_ft daxpy_blas_ptr = bli_cntx_get_l1v_ker_dt(BLIS_DOUBLE, BLIS_AXPYV_KER, cntx);

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --ipiv;

    /* Function Body */
    /* KV is the number of superdiagonals in the factor U, allowing for fill-in. */
    kv = *ku + *kl;

    /* Test the input parameters. */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*kl < 0)
    {
        *info = -3;
    }
    else if(*ku < 0)
    {
        *info = -4;
    }
    else if(*ldab < *kl + kv + 1)
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGBTF2", &i__1, (ftnlen)6);
        return;
    }

    /* Gaussian elimination with partial pivoting */
    /* Set fill-in elements in columns KU+2 to KV to zero. */
    i__1 = fla_min(kv, *n);
    for(j = *ku + 2; j <= i__1; ++j)
    {
        i__2 = *kl;
        for(i__ = kv - j + 2; i__ <= i__2; ++i__)
        {
            ab[i__ + j * ab_dim1] = 0.;
        }
    }

    /* JU is the index of the last column affected by the current stage */
    /* of the factorization. */
    ju = 1;
    i__1 = fla_min(*m, *n);
    for(j = 1; j <= i__1; ++j)
    {
        /* Set fill-in elements in column J+KV to zero. */
        if(j + kv <= *n)
        {
            i__2 = *kl;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                ab[i__ + (j + kv) * ab_dim1] = 0.;
            }
        }

        /* Find pivot and test for singularity. KM is the number of */
        /* subdiagonal elements in the current column. */
        /* Computing MIN */
        i__2 = *kl;
        i__3 = *m - j; // , expr subst
        km = fla_min(i__2, i__3);
        i__2 = km + 1;

        /* idamax_blas_ptr finds the index of the first element having maximum absolute value */
        idamax_blas_ptr((dim_t)i__2, &ab[kv + 1 + j * ab_dim1], (dim_t)c__1, &index, cntx);

        jp = (integer)index + 1;
        ipiv[j] = jp + j - 1;
        if(ab[kv + jp + j * ab_dim1] != 0.)
        {
            /* Computing MAX */
            /* Computing MIN */
            i__4 = j + *ku + jp - 1;
            i__2 = ju;
            i__3 = fla_min(i__4, *n); // , expr subst
            ju = fla_max(i__2, i__3);
            /* Apply interchange to columns J to JU. */
            if(jp != 1)
            {
                i__2 = ju - j + 1;
                i__3 = *ldab - 1;
                i__4 = *ldab - 1;

                /* dswap_blas_ptr swaps two vectors using AOCL-BLAS */
                dswap_blas_ptr((dim_t)i__2, &ab[kv + jp + j * ab_dim1], (dim_t)i__3,
                               &ab[kv + 1 + j * ab_dim1], (dim_t)i__4, cntx);
            }

            if(km > 0)
            {
                /* Compute multipliers. */
                d__1 = 1. / ab[kv + 1 + j * ab_dim1];

                /* dscal_blas_ptr scales a vector by a constant */
                dscal_blas_ptr(BLIS_NO_CONJUGATE, (dim_t)km, &d__1, &ab[kv + 2 + j * ab_dim1],
                               (dim_t)c__1, cntx);

                /* Update trailing submatrix within the band. */
                if(ju > j)
                {
                    i__2 = ju - j;
                    i__3 = *ldab - 1;
                    i__4 = *ldab - 1;
                    x = &ab[kv + 2 + j * ab_dim1];
                    y = &ab[kv + (j + 1) * ab_dim1];
                    r = &ab[kv + 1 + (j + 1) * ab_dim1];

                    for(integer i = 0; i < i__2; i++)
                    {
                        alpha = -y[i * i__4];

                        if(alpha)
                            /* daxpy_blas_ptr performs the operation y = alpha * x + y */
                            daxpy_blas_ptr(BLIS_NO_CONJUGATE, (dim_t)km, &alpha, x, (dim_t)c__1,
                                           &r[i * i__3], (dim_t)c__1, cntx);
                    }
                }
            }
        }
        else
        {
            /* If pivot is zero, set INFO to the index of the pivot */
            /* unless a zero pivot has already been found. */
            if(*info == 0)
            {
                *info = j;
            }
        }
    }
    return;
    /* End of dgbtf2_blas_ver */
}

#endif

/* > \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked
 * version of th e algorithm. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGBTF2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, KL, KU, LDAB, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION AB( LDAB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBTF2 computes an LU factorization of a real m-by-n band matrix A */
/* > using partial pivoting with row interchanges. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
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
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* > On entry, the matrix A in band storage, in rows KL+1 to */
/* > 2*KL+KU+1;
rows 1 to KL of the array need not be set. */
/* > The j-th column of A is stored in the j-th column of the */
/* > array AB as follows: */
/* > AB(kl+ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=fla_min(m,j+kl) */
/* > */
/* > On exit, details of the factorization: U is stored as an */
/* > upper triangular band matrix with KL+KU superdiagonals in */
/* > rows 1 to KL+KU+1, and the multipliers used during the */
/* > factorization are stored in rows KL+KU+2 to 2*KL+KU+1. */
/* > See below for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (fla_min(M,N)) */
/* > The pivot indices;
for 1 <= i <= fla_min(M,N), row i of the */
/* > matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = +i, U(i,i) is exactly zero. The factorization */
/* > has been completed, but the factor U is exactly */
/* > singular, and division by zero will occur if it is used */
/* > to solve a system of equations. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleGBcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The band storage scheme is illustrated by the following example, when */
/* > M = N = 6, KL = 2, KU = 1: */
/* > */
/* > On entry: On exit: */
/* > */
/* > * * * + + + * * * u14 u25 u36 */
/* > * * + + + + * * u13 u24 u35 u46 */
/* > * a12 a23 a34 a45 a56 * u12 u23 u34 u45 u56 */
/* > a11 a22 a33 a44 a55 a66 u11 u22 u33 u44 u55 u66 */
/* > a21 a32 a43 a54 a65 * m21 m32 m43 m54 m65 * */
/* > a31 a42 a53 a64 * * m31 m42 m53 m64 * * */
/* > */
/* > Array elements marked * are not used by the routine;
elements marked */
/* > + need not be set on entry, but are required by the routine to store */
/* > elements of U, because of fill-in resulting from the row */
/* > interchanges. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dgbtf2_(integer *m, integer *n, integer *kl, integer *ku, doublereal *ab, integer *ldab,
             integer *ipiv, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgbtf2 inputs: m %" FLA_IS ", n %" FLA_IS ", kl %" FLA_IS ", ku %" FLA_IS
                      ", ldab %" FLA_IS "",
                      *m, *n, *kl, *ku, *ldab);
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    /* Local variables */
    integer i__, j, km, jp, ju, kv;
#ifndef FLA_ENABLE_AOCL_BLAS
    extern /* Subroutine */
        void
        dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
              doublereal *, integer *),
        dscal_(integer *, doublereal *, doublereal *, integer *),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
#endif
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */

#if FLA_ENABLE_AOCL_BLAS
    if(*m <= 200 && *n <= 200)
    {
        dgbtf2_aocl_blas_ver(m, n, kl, ku, ab, ldab, ipiv, info);
        return;
    }
#endif

    /* KV is the number of superdiagonals in the factor U, allowing for */
    /* fill-in. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --ipiv;
#if AOCL_FLA_PROGRESS_H
    AOCL_FLA_PROGRESS_VAR;
#endif

    /* Function Body */
    kv = *ku + *kl;
    /* Test the input parameters. */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*kl < 0)
    {
        *info = -3;
    }
    else if(*ku < 0)
    {
        *info = -4;
    }
    else if(*ldab < *kl + kv + 1)
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGBTF2", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#if AOCL_FLA_PROGRESS_H
    progress_step_count = 0;
#ifndef FLA_ENABLE_WINDOWS_BUILD
    if(!aocl_fla_progress_ptr)
        aocl_fla_progress_ptr = aocl_fla_progress;
#endif
#endif

    /* Gaussian elimination with partial pivoting */
    /* Set fill-in elements in columns KU+2 to KV to zero. */
    i__1 = fla_min(kv, *n);
    for(j = *ku + 2; j <= i__1; ++j)
    {
        i__2 = *kl;
        for(i__ = kv - j + 2; i__ <= i__2; ++i__)
        {
            ab[i__ + j * ab_dim1] = 0.;
            /* L10: */
        }
        /* L20: */
    }
    /* JU is the index of the last column affected by the current stage */
    /* of the factorization. */
    ju = 1;
    i__1 = fla_min(*m, *n);
    for(j = 1; j <= i__1; ++j)
    {
#if AOCL_FLA_PROGRESS_H
        if(aocl_fla_progress_ptr)
        {
            if(j % 32 == 0 || j == i__1)
            {
                progress_step_count = j;
                AOCL_FLA_PROGRESS_FUNC_PTR("DGBTF2", 6, &progress_step_count, &progress_thread_id,
                                           &progress_total_threads);
            }
        }
#endif

        /* Set fill-in elements in column J+KV to zero. */
        if(j + kv <= *n)
        {
            i__2 = *kl;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                ab[i__ + (j + kv) * ab_dim1] = 0.;
                /* L30: */
            }
        }
        /* Find pivot and test for singularity. KM is the number of */
        /* subdiagonal elements in the current column. */
        /* Computing MIN */
        i__2 = *kl;
        i__3 = *m - j; // , expr subst
        km = fla_min(i__2, i__3);
        i__2 = km + 1;
        jp = idamax_(&i__2, &ab[kv + 1 + j * ab_dim1], &c__1);
        ipiv[j] = jp + j - 1;
        if(ab[kv + jp + j * ab_dim1] != 0.)
        {
            /* Computing MAX */
            /* Computing MIN */
            i__4 = j + *ku + jp - 1;
            i__2 = ju;
            i__3 = fla_min(i__4, *n); // , expr subst
            ju = fla_max(i__2, i__3);
            /* Apply interchange to columns J to JU. */
            if(jp != 1)
            {
                i__2 = ju - j + 1;
                i__3 = *ldab - 1;
                i__4 = *ldab - 1;
                dswap_(&i__2, &ab[kv + jp + j * ab_dim1], &i__3, &ab[kv + 1 + j * ab_dim1], &i__4);
            }
            if(km > 0)
            {
                /* Compute multipliers. */
                d__1 = 1. / ab[kv + 1 + j * ab_dim1];
                dscal_(&km, &d__1, &ab[kv + 2 + j * ab_dim1], &c__1);
                /* Update trailing submatrix within the band. */
                if(ju > j)
                {
                    i__2 = ju - j;
                    i__3 = *ldab - 1;
                    i__4 = *ldab - 1;
                    dger_(&km, &i__2, &c_b9, &ab[kv + 2 + j * ab_dim1], &c__1,
                          &ab[kv + (j + 1) * ab_dim1], &i__3, &ab[kv + 1 + (j + 1) * ab_dim1],
                          &i__4);
                }
            }
        }
        else
        {
            /* If pivot is zero, set INFO to the index of the pivot */
            /* unless a zero pivot has already been found. */
            if(*info == 0)
            {
                *info = j;
            }
        }
        /* L40: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DGBTF2 */
}
/* dgbtf2_ */