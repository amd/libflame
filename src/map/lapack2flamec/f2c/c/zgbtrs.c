/* ../netlib/zgbtrs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
/* > \brief \b ZGBTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGBTRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbtrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbtrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbtrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, KL, KU, LDAB, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 AB( LDAB, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGBTRS solves a system of linear equations */
/* > A * X = B, A**T * X = B, or A**H * X = B */
/* > with a general band matrix A using the LU factorization computed */
/* > by ZGBTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations. */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose) */
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
/* > AB is COMPLEX*16 array, dimension (LDAB,N) */
/* > Details of the LU factorization of the band matrix A, as */
/* > computed by ZGBTRF. U is stored as an upper triangular band */
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
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* > \ingroup complex16GBcomputational */
/* ===================================================================== */
/* Subroutine */
void zgbtrs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublecomplex *ab,
             integer *ldab, integer *ipiv, doublecomplex *b, integer *ldb, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgbtrs inputs: trans %c, n %" FLA_IS ", kl %" FLA_IS ", ku %" FLA_IS
                      ", nrhs %" FLA_IS ", ldab %" FLA_IS ", ldb %" FLA_IS "",
                      *trans, *n, *kl, *ku, *nrhs, *ldab, *ldb);

    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j, l, kd, lm;
    extern logical lsame_(char *, char *, integer, integer);
    logical lnoti;
    extern /* Subroutine */
        void
        zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
               doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *),
        zgeru_(integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *,
               integer *, doublecomplex *, integer *),
        zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        ztbsv_(char *, char *, char *, integer *, integer *, doublecomplex *, integer *,
               doublecomplex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        zlacgv_(integer *, doublecomplex *, integer *);
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
        xerbla_("ZGBTRS", &i__1, (ftnlen)6);
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
                    zswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgeru_(&lm, nrhs, &z__1, &ab[kd + 1 + j * ab_dim1], &c__1, &b[j + b_dim1], ldb,
                       &b[j + 1 + b_dim1], ldb);
                /* L10: */
            }
        }
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U*X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            ztbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                   &b[i__ * b_dim1 + 1], &c__1);
            /* L20: */
        }
    }
    else if(lsame_(trans, "T", 1, 1))
    {
        /* Solve A**T * X = B. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U**T * X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            ztbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                   &b[i__ * b_dim1 + 1], &c__1);
            /* L30: */
        }
        /* Solve L**T * X = B, overwriting B with X. */
        if(lnoti)
        {
            for(j = *n - 1; j >= 1; --j)
            {
                /* Computing MIN */
                i__1 = *kl;
                i__2 = *n - j; // , expr subst
                lm = fla_min(i__1, i__2);
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("Transpose", &lm, nrhs, &z__1, &b[j + 1 + b_dim1], ldb,
                       &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1, &b[j + b_dim1], ldb);
                l = ipiv[j];
                if(l != j)
                {
                    zswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                /* L40: */
            }
        }
    }
    else
    {
        /* Solve A**H * X = B. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U**H * X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            ztbsv_("Upper", "Conjugate transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                   &b[i__ * b_dim1 + 1], &c__1);
            /* L50: */
        }
        /* Solve L**H * X = B, overwriting B with X. */
        if(lnoti)
        {
            for(j = *n - 1; j >= 1; --j)
            {
                /* Computing MIN */
                i__1 = *kl;
                i__2 = *n - j; // , expr subst
                lm = fla_min(i__1, i__2);
                zlacgv_(nrhs, &b[j + b_dim1], ldb);
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("Conjugate transpose", &lm, nrhs, &z__1, &b[j + 1 + b_dim1], ldb,
                       &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1, &b[j + b_dim1], ldb);
                zlacgv_(nrhs, &b[j + b_dim1], ldb);
                l = ipiv[j];
                if(l != j)
                {
                    zswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                /* L60: */
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGBTRS */
}
/* zgbtrs_ */
