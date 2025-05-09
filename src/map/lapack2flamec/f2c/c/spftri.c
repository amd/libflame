/* ../netlib/spftri.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b11 = 1.f;
/* > \brief \b SPFTRI */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPFTRI + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spftri.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spftri.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spftri.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPFTRI( TRANSR, UPLO, N, A, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANSR, UPLO */
/* INTEGER INFO, N */
/* .. Array Arguments .. */
/* REAL A( 0: * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPFTRI computes the inverse of a real (symmetric) positive definite */
/* > matrix A using the Cholesky factorization A = U**T*U or A = L*L**T */
/* > computed by SPFTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANSR */
/* > \verbatim */
/* > TRANSR is CHARACTER*1 */
/* > = 'N': The Normal TRANSR of RFP A is stored;
 */
/* > = 'T': The Transpose TRANSR of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension ( N*(N+1)/2 ) */
/* > On entry, the symmetric matrix A in RFP format. RFP format is */
/* > described by TRANSR, UPLO, and N as follows: If TRANSR = 'N' */
/* > then RFP A is (0:N,0:k-1) when N is even;
k=N/2. RFP A is */
/* > (0:N-1,0:k) when N is odd;
k=N/2. IF TRANSR = 'T' then RFP is */
/* > the transpose of RFP A as defined when */
/* > TRANSR = 'N'. The contents of RFP A are defined by UPLO as */
/* > follows: If UPLO = 'U' the RFP A contains the nt elements of */
/* > upper packed A. If UPLO = 'L' the RFP A contains the elements */
/* > of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR = */
/* > 'T'. When TRANSR is 'N' the LDA is N+1 when N is even and N */
/* > is odd. See the Note below for more details. */
/* > */
/* > On exit, the symmetric inverse of the original matrix, in the */
/* > same storage format. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the (i,i) element of the factor U or L is */
/* > zero, and the inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > We first consider Rectangular Full Packed (RFP) Format when N is */
/* > even. We give an example where N = 6. */
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
/* > the transpose of the first three columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:2,0:2) consists of */
/* > the transpose of the last three columns of AP lower. */
/* > This covers the case N even and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > 03 04 05 33 43 53 */
/* > 13 14 15 00 44 54 */
/* > 23 24 25 10 11 55 */
/* > 33 34 35 20 21 22 */
/* > 00 44 45 30 31 32 */
/* > 01 11 55 40 41 42 */
/* > 02 12 22 50 51 52 */
/* > */
/* > Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* > RFP A RFP A */
/* > */
/* > 03 13 23 33 00 01 02 33 00 10 20 30 40 50 */
/* > 04 14 24 34 44 11 12 43 44 11 21 31 41 51 */
/* > 05 15 25 35 45 55 22 53 54 55 22 32 42 52 */
/* > */
/* > */
/* > We then consider Rectangular Full Packed (RFP) Format when N is */
/* > odd. We give an example where N = 5. */
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
/* > the transpose of the first two columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:1,1:2) consists of */
/* > the transpose of the last two columns of AP lower. */
/* > This covers the case N odd and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > 02 03 04 00 33 43 */
/* > 12 13 14 10 11 44 */
/* > 22 23 24 20 21 22 */
/* > 00 33 34 30 31 32 */
/* > 01 11 44 40 41 42 */
/* > */
/* > Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > RFP A RFP A */
/* > */
/* > 02 12 22 00 01 00 10 20 30 40 50 */
/* > 03 13 23 33 11 33 11 21 31 41 51 */
/* > 04 14 24 34 44 43 44 22 32 42 52 */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void spftri_(char *transr, char *uplo, integer *n, real *a, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spftri inputs: transr %c ,uplo %c ,n %" FLA_IS "", *transr, *uplo, *n);
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    integer k, n1, n2;
    logical normaltransr;
    extern logical lsame_(char *, char *, integer, integer);
    logical lower;
    extern /* Subroutine */
        void
        strmm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *,
               real *, integer *),
        ssyrk_(char *, char *, integer *, integer *, real *, real *, integer *, real *, real *,
               integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical nisodd;
    extern /* Subroutine */
        void
        slauum_(char *, integer *, real *, integer *, integer *),
        stftri_(char *, char *, char *, integer *, real *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
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
    *info = 0;
    normaltransr = lsame_(transr, "N", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    if(!normaltransr && !lsame_(transr, "T", 1, 1))
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
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SPFTRI", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Invert the triangular Cholesky factor U or L. */
    stftri_(transr, uplo, "N", n, a, info);
    if(*info > 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* If N is odd, set NISODD = .TRUE. */
    /* If N is even, set K = N/2 and NISODD = .FALSE. */
    if(*n % 2 == 0)
    {
        k = *n / 2;
        nisodd = FALSE_;
    }
    else
    {
        nisodd = TRUE_;
    }
    /* Set N1 and N2 depending on LOWER */
    if(lower)
    {
        n2 = *n / 2;
        n1 = *n - n2;
    }
    else
    {
        n1 = *n / 2;
        n2 = *n - n1;
    }
    /* Start execution of triangular matrix multiply: inv(U)*inv(U)^C or */
    /* inv(L)^C*inv(L). There are eight cases. */
    if(nisodd)
    {
        /* N is odd */
        if(normaltransr)
        {
            /* N is odd and TRANSR = 'N' */
            if(lower)
            {
                /* SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:N1-1) ) */
                /* T1 -> a(0,0), T2 -> a(0,1), S -> a(N1,0) */
                /* T1 -> a(0), T2 -> a(n), S -> a(N1) */
                slauum_("L", &n1, a, n, info);
                ssyrk_("L", "T", &n1, &n2, &c_b11, &a[n1], n, &c_b11, a, n);
                strmm_("L", "U", "N", "N", &n2, &n1, &c_b11, &a[*n], n, &a[n1], n);
                slauum_("U", &n2, &a[*n], n, info);
            }
            else
            {
                /* SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:N2-1) */
                /* T1 -> a(N1+1,0), T2 -> a(N1,0), S -> a(0,0) */
                /* T1 -> a(N2), T2 -> a(N1), S -> a(0) */
                slauum_("L", &n1, &a[n2], n, info);
                ssyrk_("L", "N", &n1, &n2, &c_b11, a, n, &c_b11, &a[n2], n);
                strmm_("R", "U", "T", "N", &n1, &n2, &c_b11, &a[n1], n, a, n);
                slauum_("U", &n2, &a[n1], n, info);
            }
        }
        else
        {
            /* N is odd and TRANSR = 'T' */
            if(lower)
            {
                /* SRPA for LOWER, TRANSPOSE, and N is odd */
                /* T1 -> a(0), T2 -> a(1), S -> a(0+N1*N1) */
                slauum_("U", &n1, a, &n1, info);
                ssyrk_("U", "N", &n1, &n2, &c_b11, &a[n1 * n1], &n1, &c_b11, a, &n1);
                strmm_("R", "L", "N", "N", &n1, &n2, &c_b11, &a[1], &n1, &a[n1 * n1], &n1);
                slauum_("L", &n2, &a[1], &n1, info);
            }
            else
            {
                /* SRPA for UPPER, TRANSPOSE, and N is odd */
                /* T1 -> a(0+N2*N2), T2 -> a(0+N1*N2), S -> a(0) */
                slauum_("U", &n1, &a[n2 * n2], &n2, info);
                ssyrk_("U", "T", &n1, &n2, &c_b11, a, &n2, &c_b11, &a[n2 * n2], &n2);
                strmm_("L", "L", "T", "N", &n2, &n1, &c_b11, &a[n1 * n2], &n2, a, &n2);
                slauum_("L", &n2, &a[n1 * n2], &n2, info);
            }
        }
    }
    else
    {
        /* N is even */
        if(normaltransr)
        {
            /* N is even and TRANSR = 'N' */
            if(lower)
            {
                /* SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
                /* T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
                /* T1 -> a(1), T2 -> a(0), S -> a(k+1) */
                i__1 = *n + 1;
                slauum_("L", &k, &a[1], &i__1, info);
                i__1 = *n + 1;
                i__2 = *n + 1;
                ssyrk_("L", "T", &k, &k, &c_b11, &a[k + 1], &i__1, &c_b11, &a[1], &i__2);
                i__1 = *n + 1;
                i__2 = *n + 1;
                strmm_("L", "U", "N", "N", &k, &k, &c_b11, a, &i__1, &a[k + 1], &i__2);
                i__1 = *n + 1;
                slauum_("U", &k, a, &i__1, info);
            }
            else
            {
                /* SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
                /* T1 -> a(k+1,0) , T2 -> a(k,0), S -> a(0,0) */
                /* T1 -> a(k+1), T2 -> a(k), S -> a(0) */
                i__1 = *n + 1;
                slauum_("L", &k, &a[k + 1], &i__1, info);
                i__1 = *n + 1;
                i__2 = *n + 1;
                ssyrk_("L", "N", &k, &k, &c_b11, a, &i__1, &c_b11, &a[k + 1], &i__2);
                i__1 = *n + 1;
                i__2 = *n + 1;
                strmm_("R", "U", "T", "N", &k, &k, &c_b11, &a[k], &i__1, a, &i__2);
                i__1 = *n + 1;
                slauum_("U", &k, &a[k], &i__1, info);
            }
        }
        else
        {
            /* N is even and TRANSR = 'T' */
            if(lower)
            {
                /* SRPA for LOWER, TRANSPOSE, and N is even (see paper) */
                /* T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1), */
                /* T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1));
                lda=k */
                slauum_("U", &k, &a[k], &k, info);
                ssyrk_("U", "N", &k, &k, &c_b11, &a[k * (k + 1)], &k, &c_b11, &a[k], &k);
                strmm_("R", "L", "N", "N", &k, &k, &c_b11, a, &k, &a[k * (k + 1)], &k);
                slauum_("L", &k, a, &k, info);
            }
            else
            {
                /* SRPA for UPPER, TRANSPOSE, and N is even (see paper) */
                /* T1 -> B(0,k+1), T2 -> B(0,k), S -> B(0,0), */
                /* T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0));
                lda=k */
                slauum_("U", &k, &a[k * (k + 1)], &k, info);
                ssyrk_("U", "T", &k, &k, &c_b11, a, &k, &c_b11, &a[k * (k + 1)], &k);
                strmm_("L", "L", "T", "N", &k, &k, &c_b11, &a[k * k], &k, a, &k);
                slauum_("L", &k, &a[k * k], &k, info);
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPFTRI */
}
/* spftri_ */
