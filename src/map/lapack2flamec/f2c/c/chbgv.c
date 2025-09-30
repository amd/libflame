/* ../netlib/chbgv.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CHBGST */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHBGV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbgv.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbgv.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbgv.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, */
/* LDZ, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KA, KB, LDAB, LDBB, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX AB( LDAB, * ), BB( LDBB, * ), WORK( * ), */
/* $ Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBGV computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a scomplex generalized Hermitian-definite banded eigenproblem, of */
/* > the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian */
/* > and banded, and B is also positive definite. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBZ */
/* > \verbatim */
/* > JOBZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only;
 */
/* > = 'V': Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangles of A and B are stored;
 */
/* > = 'L': Lower triangles of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KA */
/* > \verbatim */
/* > KA is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* > KB is INTEGER */
/* > The number of superdiagonals of the matrix B if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB, N) */
/* > On entry, the upper or lower triangle of the Hermitian band */
/* > matrix A, stored in the first ka+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for fla_max(1,j-ka)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+ka). */
/* > */
/* > On exit, the contents of AB are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KA+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BB */
/* > \verbatim */
/* > BB is COMPLEX array, dimension (LDBB, N) */
/* > On entry, the upper or lower triangle of the Hermitian band */
/* > matrix B, stored in the first kb+1 rows of the array. The */
/* > j-th column of B is stored in the j-th column of the array BB */
/* > as follows: */
/* > if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for fla_max(1,j-kb)<=i<=j;
 */
/* > if UPLO = 'L', BB(1+i-j,j) = B(i,j) for j<=i<=fla_min(n,j+kb). */
/* > */
/* > On exit, the factor S from the split Cholesky factorization */
/* > B = S**H*S, as returned by CPBSTF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* > LDBB is INTEGER */
/* > The leading dimension of the array BB. LDBB >= KB+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* > eigenvectors, with the i-th column of Z holding the */
/* > eigenvector associated with W(i). The eigenvectors are */
/* > normalized so that Z**H*B*Z = I. */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1, and if */
/* > JOBZ = 'V', LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, and i is: */
/* > <= N: the algorithm failed to converge: */
/* > i off-diagonal elements of an intermediate */
/* > tridiagonal form did not converge to zero;
 */
/* > > N: if INFO = N + i, for 1 <= i <= N, then CPBSTF */
/* > returned INFO = i: B is not positive definite. */
/* > The factorization of B could not be completed and */
/* > no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexOTHEReigen */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chbgv_(char *jobz, char *uplo, aocl_int_t *n, aocl_int_t *ka, aocl_int_t *kb, scomplex *ab,
            aocl_int_t *ldab, scomplex *bb, aocl_int_t *ldbb, real *w, scomplex *z__, aocl_int_t *ldz,
            scomplex *work, real *rwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z__, ldz, work, rwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ka_64 = *ka;
    aocl_int64_t kb_64 = *kb;
    aocl_int64_t ldab_64 = *ldab;
    aocl_int64_t ldbb_64 = *ldbb;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chbgv(jobz, uplo, &n_64, &ka_64, &kb_64, ab, &ldab_64, bb, &ldbb_64, w, z__,
                      &ldz_64, work, rwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_chbgv(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb,
                       scomplex *ab, aocl_int64_t *ldab, scomplex *bb, aocl_int64_t *ldbb, real *w,
                       scomplex *z__, aocl_int64_t *ldz, scomplex *work, real *rwork,
                       aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(
        buffer, 256,
        "chbgv inputs: jobz %c, uplo %c, n %lld, ka %lld, kb %lld, ldab %lld, ldbb %lld, ldz %lld",
        *jobz, *uplo, *n, *ka, *kb, *ldab, *ldbb, *ldz);
#else
    snprintf(buffer, 256,
             "chbgv inputs: jobz %c, uplo %c, n %d, ka %d, kb %d, ldab %d, ldbb %d, ldz %d", *jobz,
             *uplo, *n, *ka, *kb, *ldab, *ldbb, *ldz);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t ab_dim1, ab_offset, bb_dim1, bb_offset, z_dim1, z_offset, i__1;
    /* Local variables */
    aocl_int64_t inde;
    char vect[1];
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo;
    logical upper, wantz;
    aocl_int64_t indwrk;
    /* -- LAPACK driver routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    bb_dim1 = *ldbb;
    bb_offset = 1 + bb_dim1;
    bb -= bb_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(upper || lsame_(uplo, "L", 1, 1)))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*ka < 0)
    {
        *info = -4;
    }
    else if(*kb < 0 || *kb > *ka)
    {
        *info = -5;
    }
    else if(*ldab < *ka + 1)
    {
        *info = -7;
    }
    else if(*ldbb < *kb + 1)
    {
        *info = -9;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -12;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CHBGV ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Form a split Cholesky factorization of B. */
    aocl_lapack_cpbstf(uplo, n, kb, &bb[bb_offset], ldbb, info);
    if(*info != 0)
    {
        *info = *n + *info;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Transform problem to standard eigenvalue problem. */
    inde = 1;
    indwrk = inde + *n;
    aocl_lapack_chbgst(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
                       &z__[z_offset], ldz, &work[1], &rwork[indwrk], &iinfo);
    /* Reduce to tridiagonal form. */
    if(wantz)
    {
        *(unsigned char *)vect = 'U';
    }
    else
    {
        *(unsigned char *)vect = 'N';
    }
    aocl_lapack_chbtrd(vect, uplo, n, ka, &ab[ab_offset], ldab, &w[1], &rwork[inde], &z__[z_offset],
                       ldz, &work[1], &iinfo);
    /* For eigenvalues only, call SSTERF. For eigenvectors, call CSTEQR. */
    if(!wantz)
    {
        aocl_lapack_ssterf(n, &w[1], &rwork[inde], info);
    }
    else
    {
        aocl_lapack_csteqr(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[indwrk], info);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHBGV */
}
/* chbgv_ */
