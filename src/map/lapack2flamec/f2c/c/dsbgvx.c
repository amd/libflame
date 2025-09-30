/* ./dsbgvx.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static doublereal c_b25 = 1.;
static doublereal c_b27 = 0.;
/* > \brief \b DSBGVX */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSBGVX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgvx.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgvx.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgvx.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, */
/* LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, */
/* LDZ, WORK, IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE, UPLO */
/* INTEGER IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M, */
/* $ N */
/* DOUBLE PRECISION ABSTOL, VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IFAIL( * ), IWORK( * ) */
/* DOUBLE PRECISION AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ), */
/* $ W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBGVX computes selected eigenvalues, and optionally, eigenvectors */
/* > of a real generalized symmetric-definite banded eigenproblem, of */
/* > the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric */
/* > and banded, and B is also positive definite. Eigenvalues and */
/* > eigenvectors can be selected by specifying either all eigenvalues, */
/* > a range of values or a range of indices for the desired eigenvalues. */
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
/* > \param[in] RANGE */
/* > \verbatim */
/* > RANGE is CHARACTER*1 */
/* > = 'A': all eigenvalues will be found. */
/* > = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* > will be found. */
/* > = 'I': the IL-th through IU-th eigenvalues will be found. */
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
/* > AB is DOUBLE PRECISION array, dimension (LDAB, N) */
/* > On entry, the upper or lower triangle of the symmetric band */
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
/* > BB is DOUBLE PRECISION array, dimension (LDBB, N) */
/* > On entry, the upper or lower triangle of the symmetric band */
/* > matrix B, stored in the first kb+1 rows of the array. The */
/* > j-th column of B is stored in the j-th column of the array BB */
/* > as follows: */
/* > if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for fla_max(1,j-kb)<=i<=j;
 */
/* > if UPLO = 'L', BB(1+i-j,j) = B(i,j) for j<=i<=fla_min(n,j+kb). */
/* > */
/* > On exit, the factor S from the split Cholesky factorization */
/* > B = S**T*S, as returned by DPBSTF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* > LDBB is INTEGER */
/* > The leading dimension of the array BB. LDBB >= KB+1. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* > If JOBZ = 'V', the n-by-n matrix used in the reduction of */
/* > A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x, */
/* > and consequently C to tridiagonal form. */
/* > If JOBZ = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. If JOBZ = 'N', */
/* > LDQ >= 1. If JOBZ = 'V', LDQ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* > VL is DOUBLE PRECISION */
/* > */
/* > If RANGE='V', the lower bound of the interval to */
/* > be searched for eigenvalues. VL < VU. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* > VU is DOUBLE PRECISION */
/* > */
/* > If RANGE='V', the upper bound of the interval to */
/* > be searched for eigenvalues. VL < VU. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* > IL is INTEGER */
/* > */
/* > If RANGE='I', the index of the */
/* > smallest eigenvalue to be returned. */
/* > 1 <= IL <= IU <= N, if N > 0;
IL = 1 and IU = 0 if N = 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* > IU is INTEGER */
/* > */
/* > If RANGE='I', the index of the */
/* > largest eigenvalue to be returned. */
/* > 1 <= IL <= IU <= N, if N > 0;
IL = 1 and IU = 0 if N = 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* > ABSTOL is DOUBLE PRECISION */
/* > The absolute error tolerance for the eigenvalues. */
/* > An approximate eigenvalue is accepted as converged */
/* > when it is determined to lie in an interval [a,b] */
/* > of width less than or equal to */
/* > */
/* > ABSTOL + EPS * fla_max( |a|,|b| ) , */
/* > */
/* > where EPS is the machine precision. If ABSTOL is less than */
/* > or equal to zero, then EPS*|T| will be used in its place, */
/* > where |T| is the 1-norm of the tridiagonal matrix obtained */
/* > by reducing A to tridiagonal form. */
/* > */
/* > Eigenvalues will be computed most accurately when ABSTOL is */
/* > set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
/* > If this routine returns with INFO>0, indicating that some */
/* > eigenvectors did not converge, try setting ABSTOL to */
/* > 2*DLAMCH('S'). */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The total number of eigenvalues found. 0 <= M <= N. */
/* > If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is DOUBLE PRECISION array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* > eigenvectors, with the i-th column of Z holding the */
/* > eigenvector associated with W(i). The eigenvectors are */
/* > normalized so Z**T*B*Z = I. */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1, and if */
/* > JOBZ = 'V', LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (7*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* > IFAIL is INTEGER array, dimension (M) */
/* > If JOBZ = 'V', then if INFO = 0, the first M elements of */
/* > IFAIL are zero. If INFO > 0, then IFAIL contains the */
/* > indices of the eigenvalues that failed to converge. */
/* > If JOBZ = 'N', then IFAIL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > <= N: if INFO = i, then i eigenvectors failed to converge. */
/* > Their indices are stored in IFAIL. */
/* > > N: DPBSTF returned an error code;
i.e., */
/* > if INFO = N + i, for 1 <= i <= N, then the leading */
/* > principal minor of order i of B is not positive. */
/* > The factorization of B could not be completed and */
/* > no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup hbgvx */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dsbgvx_(char *jobz, char *range, char *uplo, aocl_int_t *n, aocl_int_t *ka, aocl_int_t *kb,
             doublereal *ab, aocl_int_t *ldab, doublereal *bb, aocl_int_t *ldbb, doublereal *q,
             aocl_int_t *ldq, doublereal *vl, doublereal *vu, aocl_int_t *il, aocl_int_t *iu,
             doublereal *abstol, aocl_int_t *m, doublereal *w, doublereal *z__, aocl_int_t *ldz,
             doublereal *work, aocl_int_t *iwork, aocl_int_t *ifail, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu,
                       abstol, m, w, z__, ldz, work, iwork, ifail, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ka_64 = *ka;
    aocl_int64_t kb_64 = *kb;
    aocl_int64_t ldab_64 = *ldab;
    aocl_int64_t ldbb_64 = *ldbb;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t il_64 = *il;
    aocl_int64_t iu_64 = *iu;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsbgvx(jobz, range, uplo, &n_64, &ka_64, &kb_64, ab, &ldab_64, bb, &ldbb_64, q,
                       &ldq_64, vl, vu, &il_64, &iu_64, abstol, &m_64, w, z__, &ldz_64, work, iwork,
                       ifail, &info_64);

    *m = (aocl_int_t)m_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_dsbgvx(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *ka,
                        aocl_int64_t *kb, doublereal *ab, aocl_int64_t *ldab, doublereal *bb,
                        aocl_int64_t *ldbb, doublereal *q, aocl_int64_t *ldq, doublereal *vl,
                        doublereal *vu, aocl_int64_t *il, aocl_int64_t *iu, doublereal *abstol,
                        aocl_int64_t *m, doublereal *w, doublereal *z__, aocl_int64_t *ldz,
                        doublereal *work, aocl_int_t *iwork, aocl_int_t *ifail, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dsbgvx inputs: jobz %c, range %c, uplo %c, n %" FLA_IS ", ka %" FLA_IS
                      ", kb %" FLA_IS ", ldab %" FLA_IS ", ldbb %" FLA_IS ", ldq %" FLA_IS
                      ", il %" FLA_IS ", iu %" FLA_IS ", ldz %" FLA_IS "",
                      *jobz, *range, *uplo, *n, *ka, *kb, *ldab, *ldbb, *ldq, *il, *iu, *ldz);
    /* System generated locals */
    aocl_int64_t ab_dim1, ab_offset, bb_dim1, bb_offset, q_dim1, q_offset, z_dim1, z_offset, i__1,
        i__2;
    /* Local variables */
    aocl_int64_t i__, j, jj;
    doublereal tmp1;
    aocl_int64_t indd, inde;
    char vect[1];
    logical test;
    aocl_int64_t itmp1, indee;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo;
    char order[1];
    logical upper, wantz, alleig, indeig, valeig;
    aocl_int64_t indisp;
    aocl_int64_t indiwo;
    aocl_int64_t indwrk;
    aocl_int64_t nsplit;
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    bb_dim1 = *ldbb;
    bb_offset = 1 + bb_dim1;
    bb -= bb_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    alleig = lsame_(range, "A", 1, 1);
    valeig = lsame_(range, "V", 1, 1);
    indeig = lsame_(range, "I", 1, 1);
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(alleig || valeig || indeig))
    {
        *info = -2;
    }
    else if(!(upper || lsame_(uplo, "L", 1, 1)))
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*ka < 0)
    {
        *info = -5;
    }
    else if(*kb < 0 || *kb > *ka)
    {
        *info = -6;
    }
    else if(*ldab < *ka + 1)
    {
        *info = -8;
    }
    else if(*ldbb < *kb + 1)
    {
        *info = -10;
    }
    else if(*ldq < 1 || wantz && *ldq < *n)
    {
        *info = -12;
    }
    else
    {
        if(valeig)
        {
            if(*n > 0 && *vu <= *vl)
            {
                *info = -14;
            }
        }
        else if(indeig)
        {
            if(*il < 1 || *il > fla_max(1, *n))
            {
                *info = -15;
            }
            else if(*iu < fla_min(*n, *il) || *iu > *n)
            {
                *info = -16;
            }
        }
    }
    if(*info == 0)
    {
        if(*ldz < 1 || wantz && *ldz < *n)
        {
            *info = -21;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DSBGVX", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *m = 0;
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Form a split Cholesky factorization of B. */
    aocl_lapack_dpbstf(uplo, n, kb, &bb[bb_offset], ldbb, info);
    if(*info != 0)
    {
        *info = *n + *info;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Transform problem to standard eigenvalue problem. */
    aocl_lapack_dsbgst(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
                       &q[q_offset], ldq, &work[1], &iinfo);
    /* Reduce symmetric band matrix to tridiagonal form. */
    indd = 1;
    inde = indd + *n;
    indwrk = inde + *n;
    if(wantz)
    {
        *(unsigned char *)vect = 'U';
    }
    else
    {
        *(unsigned char *)vect = 'N';
    }
    aocl_lapack_dsbtrd(vect, uplo, n, ka, &ab[ab_offset], ldab, &work[indd], &work[inde],
                       &q[q_offset], ldq, &work[indwrk], &iinfo);
    /* If all eigenvalues are desired and ABSTOL is less than or equal */
    /* to zero, then call DSTERF or SSTEQR. If this fails for some */
    /* eigenvalue, then try DSTEBZ. */
    test = FALSE_;
    if(indeig)
    {
        if(*il == 1 && *iu == *n)
        {
            test = TRUE_;
        }
    }
    if((alleig || test) && *abstol <= 0.)
    {
        aocl_blas_dcopy(n, &work[indd], &c__1, &w[1], &c__1);
        indee = indwrk + (*n << 1);
        i__1 = *n - 1;
        aocl_blas_dcopy(&i__1, &work[inde], &c__1, &work[indee], &c__1);
        if(!wantz)
        {
            aocl_lapack_dsterf(n, &w[1], &work[indee], info);
        }
        else
        {
            aocl_lapack_dlacpy("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz);
            aocl_lapack_dsteqr(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[indwrk],
                               info);
            if(*info == 0)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    ifail[i__] = 0;
                    /* L10: */
                }
            }
        }
        if(*info == 0)
        {
            *m = *n;
            goto L30;
        }
        *info = 0;
    }
    /* Otherwise, call DSTEBZ and, if eigenvectors are desired, */
    /* call DSTEIN. */
    if(wantz)
    {
        *(unsigned char *)order = 'B';
    }
    else
    {
        *(unsigned char *)order = 'E';
    }
    indisp = *n + 1;
    indiwo = indisp + *n;
    aocl_lapack_dstebz(range, order, n, vl, vu, il, iu, abstol, &work[indd], &work[inde], m,
                       &nsplit, &w[1], &iwork[1], &iwork[indisp], &work[indwrk], &iwork[indiwo],
                       info);
    if(wantz)
    {
        aocl_lapack_dstein(n, &work[indd], &work[inde], m, &w[1], &iwork[1], &iwork[indisp],
                           &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &ifail[1], info);
        /* Apply transformation matrix used in reduction to tridiagonal */
        /* form to eigenvectors returned by DSTEIN. */
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            aocl_blas_dcopy(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
            aocl_blas_dgemv("N", n, n, &c_b25, &q[q_offset], ldq, &work[1], &c__1, &c_b27,
                            &z__[j * z_dim1 + 1], &c__1);
            /* L20: */
        }
    }
L30: /* If eigenvalues are not in order, then sort them, along with */
    /* eigenvectors. */
    if(wantz)
    {
        i__1 = *m - 1;
        for(j = 1; j <= i__1; ++j)
        {
            i__ = 0;
            tmp1 = w[j];
            i__2 = *m;
            for(jj = j + 1; jj <= i__2; ++jj)
            {
                if(w[jj] < tmp1)
                {
                    i__ = jj;
                    tmp1 = w[jj];
                }
                /* L40: */
            }
            if(i__ != 0)
            {
                itmp1 = iwork[i__];
                w[i__] = w[j];
                iwork[i__] = iwork[j];
                w[j] = tmp1;
                iwork[j] = (aocl_int_t)(itmp1);
                aocl_blas_dswap(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
                if(*info != 0)
                {
                    itmp1 = ifail[i__];
                    ifail[i__] = ifail[j];
                    ifail[j] = (aocl_int_t)(itmp1);
                }
            }
            /* L50: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DSBGVX */
}
/* dsbgvx_ */
