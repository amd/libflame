/* ../netlib/cpftrs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
/* > \brief \b CPFTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPFTRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpftrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpftrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpftrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPFTRS( TRANSR, UPLO, N, NRHS, A, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANSR, UPLO */
/* INTEGER INFO, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( 0: * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPFTRS solves a system of linear equations A*X = B with a Hermitian */
/* > positive definite matrix A using the Cholesky factorization */
/* > A = U**H*U or A = L*L**H computed by CPFTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANSR */
/* > \verbatim */
/* > TRANSR is CHARACTER*1 */
/* > = 'N': The Normal TRANSR of RFP A is stored;
 */
/* > = 'C': The Conjugate-transpose TRANSR of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of RFP A is stored;
 */
/* > = 'L': Lower triangle of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension ( N*(N+1)/2 );
 */
/* > The triangular factor U or L from the Cholesky factorization */
/* > of RFP A = U**H*U or RFP A = L*L**H, as computed by CPFTRF. */
/* > See note below for more details about RFP A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
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
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > We first consider Standard Packed Format when N is even. */
/* > We give an example where N = 6. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 05 00 */
/* > 11 12 13 14 15 10 11 */
/* > 22 23 24 25 20 21 22 */
/* > 33 34 35 30 31 32 33 */
/* > 44 45 40 41 42 43 44 */
/* > 55 50 51 52 53 54 55 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(4:6,0:2) consists of */
/* > conjugate-transpose of the first three columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:2,0:2) consists of */
/* > conjugate-transpose of the last three columns of AP lower. */
/* > To denote conjugate we place -- above the element. This covers the */
/* > case N even and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- -- */
/* > 03 04 05 33 43 53 */
/* > -- -- */
/* > 13 14 15 00 44 54 */
/* > -- */
/* > 23 24 25 10 11 55 */
/* > */
/* > 33 34 35 20 21 22 */
/* > -- */
/* > 00 44 45 30 31 32 */
/* > -- -- */
/* > 01 11 55 40 41 42 */
/* > -- -- -- */
/* > 02 12 22 50 51 52 */
/* > */
/* > Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate- */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- -- -- -- -- -- -- -- -- */
/* > 03 13 23 33 00 01 02 33 00 10 20 30 40 50 */
/* > -- -- -- -- -- -- -- -- -- -- */
/* > 04 14 24 34 44 11 12 43 44 11 21 31 41 51 */
/* > -- -- -- -- -- -- -- -- -- -- */
/* > 05 15 25 35 45 55 22 53 54 55 22 32 42 52 */
/* > */
/* > */
/* > We next consider Standard Packed Format when N is odd. */
/* > We give an example where N = 5. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 00 */
/* > 11 12 13 14 10 11 */
/* > 22 23 24 20 21 22 */
/* > 33 34 30 31 32 33 */
/* > 44 40 41 42 43 44 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(3:4,0:1) consists of */
/* > conjugate-transpose of the first two columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:1,1:2) consists of */
/* > conjugate-transpose of the last two columns of AP lower. */
/* > To denote conjugate we place -- above the element. This covers the */
/* > case N odd and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- */
/* > 02 03 04 00 33 43 */
/* > -- */
/* > 12 13 14 10 11 44 */
/* > */
/* > 22 23 24 20 21 22 */
/* > -- */
/* > 00 33 34 30 31 32 */
/* > -- -- */
/* > 01 11 44 40 41 42 */
/* > */
/* > Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate- */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- -- -- -- -- -- -- -- */
/* > 02 12 22 00 01 00 10 20 30 40 50 */
/* > -- -- -- -- -- -- -- -- -- */
/* > 03 13 23 33 11 33 11 21 31 41 51 */
/* > -- -- -- -- -- -- -- -- -- */
/* > 04 14 24 34 44 43 44 22 32 42 52 */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void cpftrs_(char *transr, char *uplo, integer *n, integer *nrhs, complex *a, complex *b,
             integer *ldb, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cpftrs inputs: transr %c, uplo %c, n %lld, nrhs %lld, ldb %lld", *transr,
             *uplo, *n, *nrhs, *ldb);
#else
    snprintf(buffer, 256, "cpftrs inputs: transr %c, uplo %c, n %d, nrhs %d, ldb %d", *transr,
             *uplo, *n, *nrhs, *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer b_dim1, b_offset, i__1;
    /* Local variables */
    logical normaltransr;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        ctfsm_(char *, char *, char *, char *, char *, integer *, integer *, complex *, complex *,
               complex *, integer *);
    logical lower;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    normaltransr = lsame_(transr, "N", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    if(!normaltransr && !lsame_(transr, "C", 1, 1))
    {
        *info = -1;
    }
    else if(!lower && !lsame_(uplo, "U", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*nrhs < 0)
    {
        *info = -4;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CPFTRS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* start execution: there are two triangular solves */
    if(lower)
    {
        ctfsm_(transr, "L", uplo, "N", "N", n, nrhs, &c_b1, a, &b[b_offset], ldb);
        ctfsm_(transr, "L", uplo, "C", "N", n, nrhs, &c_b1, a, &b[b_offset], ldb);
    }
    else
    {
        ctfsm_(transr, "L", uplo, "C", "N", n, nrhs, &c_b1, a, &b[b_offset], ldb);
        ctfsm_(transr, "L", uplo, "N", "N", n, nrhs, &c_b1, a, &b[b_offset], ldb);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CPFTRS */
}
/* cpftrs_ */
