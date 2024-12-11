/* ../netlib/dgbtrs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b23 = 1.;

#if FLA_ENABLE_AOCL_BLAS

/* This function is an implementation of dgbtrs using AOCL-BLAS compute kernels */
void dgbtrs_aocl_blas_ver(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                          doublereal *ab, integer *ldab, integer *ipiv, doublereal *b, integer *ldb,
                          integer *info)
{
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;
    integer i__, j, l, kd, lm;
    logical lnoti;
    logical notran;
    doublereal alpha;
    doublereal *x, *y, *r;

    /* Make a copy of AOCL-BLAS framework context. This information is needed to query the
     * architecture specific details of compute kernel */
    cntx_t *cntx = bli_gks_query_cntx();

    /* Query names of compute kernel from AOCL-BLAS framework context */
    dswapv_ker_ft dswap_blas_ptr = bli_cntx_get_l1v_ker_dt(BLIS_DOUBLE, BLIS_SWAPV_KER, cntx);
    daxpyv_ker_ft daxpy_blas_ptr = bli_cntx_get_l1v_ker_dt(BLIS_DOUBLE, BLIS_AXPYV_KER, cntx);

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", 1, 1);
    if(!notran && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
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
    else if(*nrhs < 0)
    {
        *info = -5;
    }
    else if(*ldab < (*kl << 1) + *ku + 1)
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -10;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGBTRS", &i__1, (ftnlen)6);
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        return;
    }
    kd = *ku + *kl + 1;
    lnoti = *kl > 0;
    if(notran)
    {
        /* Solve A*X = B. */
        if(lnoti)
        {
            i__1 = *n - 1;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MIN */
                i__2 = *kl;
                i__3 = *n - j; // , expr subst
                lm = fla_min(i__2, i__3);
                i__3 = *ldb;
                l = ipiv[j];
                if(l != j)
                {
                    /* dswap_blas_ptr swaps two vectors using AOCL-BLAS */
                    dswap_blas_ptr((dim_t)*nrhs, &b[l + b_dim1], (dim_t)*ldb, &b[j + b_dim1],
                                   (dim_t)*ldb, cntx);
                }

                x = &ab[kd + 1 + j * ab_dim1];
                y = &b[j + b_dim1];
                r = &b[j + 1 + b_dim1];

                for(integer i = 0; i < *nrhs; i++)
                {
                    alpha = -y[i * i__3];

                    if(alpha)
                        /* daxpy_blas_ptr performs the operation y = alpha * x + y */
                        daxpy_blas_ptr(BLIS_NO_CONJUGATE, (dim_t)lm, &alpha, x, (dim_t)c__1,
                                       &r[i * i__3], (dim_t)c__1, cntx);
                }
            }
        }
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U*X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            dtbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                   &b[i__ * b_dim1 + 1], &c__1);
            /* L20: */
        }
    }
    else
    {
        /* Solve A**T*X = B. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U**T*X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            dtbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                   &b[i__ * b_dim1 + 1], &c__1);
        }
        /* Solve L**T*X = B, overwriting B with X. */
        if(lnoti)
        {
            for(j = *n - 1; j >= 1; --j)
            {
                /* Computing MIN */
                i__1 = *kl;
                i__2 = *n - j; // , expr subst
                lm = fla_min(i__1, i__2);
                dgemv_("Transpose", &lm, nrhs, &c_b7, &b[j + 1 + b_dim1], ldb,
                       &ab[kd + 1 + j * ab_dim1], &c__1, &c_b23, &b[j + b_dim1], ldb);
                l = ipiv[j];
                if(l != j)
                {
                    dswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
            }
        }
    }
    return;
}

#endif

/* > \brief \b DGBTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGBTRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, KL, KU, LDAB, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION AB( LDAB, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBTRS solves a system of linear equations */
/* > A * X = B or A**T * X = B */
/* > with a general band matrix A using the LU factorization computed */
/* > by DGBTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations. */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T* X = B (Transpose) */
/* > = 'C': A**T* X = B (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
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
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* > Details of the LU factorization of the band matrix A, as */
/* > computed by DGBTRF. U is stored as an upper triangular band */
/* > matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and */
/* > the multipliers used during the factorization are stored in */
/* > rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices;
for 1 <= i <= N, row i of the matrix was */
/* > interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
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
/* > \date November 2011 */
/* > \ingroup doubleGBcomputational */
/* ===================================================================== */
/* Subroutine */
void dgbtrs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublereal *ab,
             integer *ldab, integer *ipiv, doublereal *b, integer *ldb, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgbtrs inputs: trans %c, n %" FLA_IS ", kl %" FLA_IS ", ku %" FLA_IS
                      ", nrhs %" FLA_IS ", ldab %" FLA_IS ", ldb %" FLA_IS "",
                      *trans, *n, *kl, *ku, *nrhs, *ldab, *ldb);
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j, l, kd, lm;
#ifndef FLA_ENABLE_AOCL_BLAS
    extern /* Subroutine */
        void
        dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
              doublereal *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *),
        dtbsv_(char *, char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
               integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
#endif

    logical lnoti;
    logical notran;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    if(*n <= 200)
    {
        dgbtrs_aocl_blas_ver(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
        return;
    }
#endif

    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", 1, 1);
    if(!notran && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
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
    else if(*nrhs < 0)
    {
        *info = -5;
    }
    else if(*ldab < (*kl << 1) + *ku + 1)
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -10;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGBTRS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    kd = *ku + *kl + 1;
    lnoti = *kl > 0;
    if(notran)
    {
        /* Solve A*X = B. */
        /* Solve L*X = B, overwriting B with X. */
        /* L is represented as a product of permutations and unit lower */
        /* triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), */
        /* where each transformation L(i) is a rank-one modification of */
        /* the identity matrix. */
        if(lnoti)
        {
            i__1 = *n - 1;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MIN */
                i__2 = *kl;
                i__3 = *n - j; // , expr subst
                lm = fla_min(i__2, i__3);
                l = ipiv[j];
                if(l != j)
                {
                    dswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                dger_(&lm, nrhs, &c_b7, &ab[kd + 1 + j * ab_dim1], &c__1, &b[j + b_dim1], ldb,
                      &b[j + 1 + b_dim1], ldb);
                /* L10: */
            }
        }
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U*X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            dtbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                   &b[i__ * b_dim1 + 1], &c__1);
            /* L20: */
        }
    }
    else
    {
        /* Solve A**T*X = B. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U**T*X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            dtbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                   &b[i__ * b_dim1 + 1], &c__1);
            /* L30: */
        }
        /* Solve L**T*X = B, overwriting B with X. */
        if(lnoti)
        {
            for(j = *n - 1; j >= 1; --j)
            {
                /* Computing MIN */
                i__1 = *kl;
                i__2 = *n - j; // , expr subst
                lm = fla_min(i__1, i__2);
                dgemv_("Transpose", &lm, nrhs, &c_b7, &b[j + 1 + b_dim1], ldb,
                       &ab[kd + 1 + j * ab_dim1], &c__1, &c_b23, &b[j + b_dim1], ldb);
                l = ipiv[j];
                if(l != j)
                {
                    dswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                /* L40: */
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DGBTRS */
}
/* dgbtrs_ */
