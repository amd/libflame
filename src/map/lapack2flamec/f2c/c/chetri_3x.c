/* ../netlib/v3.9.0/chetri_3x.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
static complex c_b2 = {0.f, 0.f};
/* > \brief \b CHETRI_3X */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHETRI_3X + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri_
 * 3x.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri_
 * 3x.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri_
 * 3x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHETRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), E( * ), WORK( N+NB+1, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > CHETRI_3X computes the inverse of a complex Hermitian indefinite */
/* > matrix A using the factorization computed by CHETRF_RK or CHETRF_BK: */
/* > */
/* > A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is Hermitian and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are */
/* > stored as an upper or lower triangular matrix. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, diagonal of the block diagonal matrix D and */
/* > factors U or L as computed by CHETRF_RK and CHETRF_BK: */
/* > a) ONLY diagonal elements of the Hermitian block diagonal */
/* > matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
 */
/* > (superdiagonal (or subdiagonal) elements of D */
/* > should be provided on entry in array E), and */
/* > b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* > If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* > On exit, if INFO = 0, the Hermitian inverse of the original */
/* > matrix. */
/* > If UPLO = 'U': the upper triangular part of the inverse */
/* > is formed and the part of A below the diagonal is not */
/* > referenced;
 */
/* > If UPLO = 'L': the lower triangular part of the inverse */
/* > is formed and the part of A above the diagonal is not */
/* > referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is COMPLEX array, dimension (N) */
/* > On entry, contains the superdiagonal (or subdiagonal) */
/* > elements of the Hermitian block diagonal matrix D */
/* > with 1-by-1 or 2-by-2 diagonal blocks, where */
/* > If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) not referenced;
 */
/* > If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) not referenced. */
/* > */
/* > NOTE: For 1-by-1 diagonal block D(k), where */
/* > 1 <= k <= N, the element E(k) is not referenced in both */
/* > UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by CHETRF_RK or CHETRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N+NB+1,NB+3). */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > Block size. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) = 0;
the matrix is singular and its */
/* > inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup complexHEcomputational */
/* > \par Contributors: */
/* ================== */
/* > \verbatim */
/* > */
/* > June 2017, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
void chetri_3x_(char *uplo, integer *n, complex *a, integer *lda, complex *e, integer *ipiv,
                complex *work, integer *nb, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "chetri_3x inputs: uplo %c, n %lld, lda %lld, nb %lld", *uplo, *n, *lda,
             *nb);
#else
    snprintf(buffer, 256, "chetri_3x inputs: uplo %c, n %d, lda %d, nb %d", *uplo, *n, *lda, *nb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    double c_abs(complex *);
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);
    /* Local variables */
    extern /* Subroutine */
        void
        cheswapr_(char *, integer *, complex *, integer *, integer *, integer *);
    complex d__;
    integer i__, j, k;
    real t, ak;
    integer u11, ip, nnb, cut;
    real akp1;
    integer invd;
    complex akkp1;
    extern /* Subroutine */
        void
        cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        ctrmm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *,
               integer *, complex *, integer *);
    logical upper;
    complex u01_i_j__, u11_i_j__;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    integer icount;
    extern /* Subroutine */
        void
        ctrtri_(char *, char *, integer *, complex *, integer *, integer *);
    complex u01_ip1_j__, u11_ip1_j__;
    /* -- LAPACK computational routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;
    work_dim1 = *n + *nb + 1;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    /* Quick return if possible */
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHETRI_3X", &i__1, (ftnlen)9);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Workspace got Non-diag elements of D */
    i__1 = *n;
    for(k = 1; k <= i__1; ++k)
    {
        i__2 = k + work_dim1;
        i__3 = k;
        work[i__2].r = e[i__3].r;
        work[i__2].i = e[i__3].i; // , expr subst
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if(upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for(*info = *n; *info >= 1; --(*info))
        {
            i__1 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__1].r == 0.f && a[i__1].i == 0.f))
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        i__1 = *n;
        for(*info = 1; *info <= i__1; ++(*info))
        {
            i__2 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__2].r == 0.f && a[i__2].i == 0.f))
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
        }
    }
    *info = 0;
    /* Splitting Workspace */
    /* U01 is a block ( N, NB+1 ) */
    /* The first element of U01 is in WORK( 1, 1 ) */
    /* U11 is a block ( NB+1, NB+1 ) */
    /* The first element of U11 is in WORK( N+1, 1 ) */
    u11 = *n;
    /* INVD is a block ( N, 2 ) */
    /* The first element of INVD is in WORK( 1, INVD ) */
    invd = *nb + 2;
    if(upper)
    {
        /* Begin Upper */
        /* invA = P * inv(U**H) * inv(D) * inv(U) * P**T. */
        ctrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D) * inv(U) */
        k = 1;
        while(k <= *n)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                i__2 = k + k * a_dim1;
                r__1 = 1.f / a[i__2].r;
                work[i__1].r = r__1;
                work[i__1].i = 0.f; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].r = 0.f;
                work[i__1].i = 0.f; // , expr subst
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                t = c_abs(&work[k + 1 + work_dim1]);
                i__1 = k + k * a_dim1;
                ak = a[i__1].r / t;
                i__1 = k + 1 + (k + 1) * a_dim1;
                akp1 = a[i__1].r / t;
                i__1 = k + 1 + work_dim1;
                q__1.r = work[i__1].r / t;
                q__1.i = work[i__1].i / t; // , expr subst
                akkp1.r = q__1.r;
                akkp1.i = q__1.i; // , expr subst
                r__1 = ak * akp1;
                q__2.r = r__1 - 1.f;
                q__2.i = -0.f; // , expr subst
                q__1.r = t * q__2.r;
                q__1.i = t * q__2.i; // , expr subst
                d__.r = q__1.r;
                d__.i = q__1.i; // , expr subst
                i__1 = k + invd * work_dim1;
                q__2.r = akp1;
                q__2.i = 0.f; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + 1 + (invd + 1) * work_dim1;
                q__2.r = ak;
                q__2.i = 0.f; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                q__2.r = -akkp1.r;
                q__2.i = -akkp1.i; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + 1 + invd * work_dim1;
                r_cnjg(&q__1, &work[k + (invd + 1) * work_dim1]);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                ++k;
            }
            ++k;
        }
        /* inv(U**H) = (inv(U))**H */
        /* inv(U**H) * inv(D) * inv(U) */
        cut = *n;
        while(cut > 0)
        {
            nnb = *nb;
            if(cut <= nnb)
            {
                nnb = cut;
            }
            else
            {
                icount = 0;
                /* count negative elements, */
                i__1 = cut;
                for(i__ = cut + 1 - nnb; i__ <= i__1; ++i__)
                {
                    if(ipiv[i__] < 0)
                    {
                        ++icount;
                    }
                }
                /* need a even number for a clear cut */
                if(icount % 2 == 1)
                {
                    ++nnb;
                }
            }
            cut -= nnb;
            /* U01 Block */
            i__1 = cut;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * work_dim1;
                    i__4 = i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
                }
            }
            /* U11 Block */
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = u11 + i__ + i__ * work_dim1;
                work[i__2].r = 1.f;
                work[i__2].i = 0.f; // , expr subst
                i__2 = i__ - 1;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].r = 0.f;
                    work[i__3].i = 0.f; // , expr subst
                }
                i__2 = nnb;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
                }
            }
            /* invD * U01 */
            i__ = 1;
            while(i__ <= cut)
            {
                if(ipiv[i__] > 0)
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        i__3 = i__ + invd * work_dim1;
                        i__4 = i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[i__4].i
                                 + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                }
                else
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        u01_i_j__.r = work[i__2].r;
                        u01_i_j__.i = work[i__2].i; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        u01_ip1_j__.r = work[i__2].r;
                        u01_ip1_j__.i = work[i__2].i; // , expr subst
                        i__2 = i__ + j * work_dim1;
                        i__3 = i__ + invd * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i
                                 + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = i__ + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i
                                 + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        i__3 = i__ + 1 + invd * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i
                                 + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = i__ + 1 + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i
                                 + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                    ++i__;
                }
                ++i__;
            }
            /* invD1 * U11 */
            i__ = 1;
            while(i__ <= nnb)
            {
                if(ipiv[cut + i__] > 0)
                {
                    i__1 = nnb;
                    for(j = i__; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[i__4].i
                                 + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                }
                else
                {
                    i__1 = nnb;
                    for(j = i__; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        u11_i_j__.r = work[i__2].r;
                        u11_i_j__.i = work[i__2].i; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        u11_ip1_j__.r = work[i__2].r;
                        u11_ip1_j__.i = work[i__2].i; // , expr subst
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__2.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__2.i = work[i__3].r * work[i__4].i
                                 + work[i__3].i * work[i__4].r; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        i__6 = u11 + i__ + 1 + j * work_dim1;
                        q__3.r = work[i__5].r * work[i__6].r - work[i__5].i * work[i__6].i;
                        q__3.i = work[i__5].r * work[i__6].i
                                 + work[i__5].i * work[i__6].r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        i__3 = cut + i__ + 1 + invd * work_dim1;
                        q__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * u11_i_j__.i;
                        q__2.i = work[i__3].r * u11_i_j__.i
                                 + work[i__3].i * u11_i_j__.r; // , expr subst
                        i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i * u11_ip1_j__.i;
                        q__3.i = work[i__4].r * u11_ip1_j__.i
                                 + work[i__4].i * u11_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                    ++i__;
                }
                ++i__;
            }
            /* U11**H * invD1 * U11 -> U11 */
            i__1 = *n + *nb + 1;
            ctrmm_("L", "U", "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1], lda,
                   &work[u11 + 1 + work_dim1], &i__1);
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = i__; j <= i__2; ++j)
                {
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = u11 + i__ + j * work_dim1;
                    a[i__3].r = work[i__4].r;
                    a[i__3].i = work[i__4].i; // , expr subst
                }
            }
            /* U01**H * invD * U01 -> A( CUT+I, CUT+J ) */
            i__1 = *n + *nb + 1;
            i__2 = *n + *nb + 1;
            cgemm_("C", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 1], lda,
                   &work[work_offset], &i__1, &c_b2, &work[u11 + 1 + work_dim1], &i__2);
            /* U11 = U11**H * invD1 * U11 + U01**H * invD * U01 */
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = i__; j <= i__2; ++j)
                {
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    i__5 = u11 + i__ + j * work_dim1;
                    q__1.r = a[i__4].r + work[i__5].r;
                    q__1.i = a[i__4].i + work[i__5].i; // , expr subst
                    a[i__3].r = q__1.r;
                    a[i__3].i = q__1.i; // , expr subst
                }
            }
            /* U01 = U00**H * invD0 * U01 */
            i__1 = *n + *nb + 1;
            ctrmm_("L", uplo, "C", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, &work[work_offset],
                   &i__1);
            /* Update U01 */
            i__1 = cut;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + (cut + j) * a_dim1;
                    i__4 = i__ + j * work_dim1;
                    a[i__3].r = work[i__4].r;
                    a[i__3].i = work[i__4].i; // , expr subst
                }
            }
            /* Next Block */
        }
        /* Apply PERMUTATIONS P and P**T: */
        /* P * inv(U**H) * inv(D) * inv(U) * P**T. */
        /* Interchange rows and columns I and IPIV(I) in reverse order */
        /* from the formation order of IPIV vector for Upper case. */
        /* ( We can use a loop over IPIV with increment 1, */
        /* since the ABS value of IPIV(I) represents the row (column) */
        /* index of the interchange with row (column) i in both 1x1 */
        /* and 2x2 pivot cases, i.e. we don't need separate code branches */
        /* for 1x1 and 2x2 pivot cases ) */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            ip = (i__2 = ipiv[i__], f2c_abs(i__2));
            if(ip != i__)
            {
                if(i__ < ip)
                {
                    cheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
        }
    }
    else
    {
        /* Begin Lower */
        /* inv A = P * inv(L**H) * inv(D) * inv(L) * P**T. */
        ctrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D) * inv(L) */
        k = *n;
        while(k >= 1)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                i__2 = k + k * a_dim1;
                r__1 = 1.f / a[i__2].r;
                work[i__1].r = r__1;
                work[i__1].i = 0.f; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].r = 0.f;
                work[i__1].i = 0.f; // , expr subst
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                t = c_abs(&work[k - 1 + work_dim1]);
                i__1 = k - 1 + (k - 1) * a_dim1;
                ak = a[i__1].r / t;
                i__1 = k + k * a_dim1;
                akp1 = a[i__1].r / t;
                i__1 = k - 1 + work_dim1;
                q__1.r = work[i__1].r / t;
                q__1.i = work[i__1].i / t; // , expr subst
                akkp1.r = q__1.r;
                akkp1.i = q__1.i; // , expr subst
                r__1 = ak * akp1;
                q__2.r = r__1 - 1.f;
                q__2.i = -0.f; // , expr subst
                q__1.r = t * q__2.r;
                q__1.i = t * q__2.i; // , expr subst
                d__.r = q__1.r;
                d__.i = q__1.i; // , expr subst
                i__1 = k - 1 + invd * work_dim1;
                q__2.r = akp1;
                q__2.i = 0.f; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + invd * work_dim1;
                q__2.r = ak;
                q__2.i = 0.f; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                q__2.r = -akkp1.r;
                q__2.i = -akkp1.i; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k - 1 + (invd + 1) * work_dim1;
                r_cnjg(&q__1, &work[k + (invd + 1) * work_dim1]);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                --k;
            }
            --k;
        }
        /* inv(L**H) = (inv(L))**H */
        /* inv(L**H) * inv(D) * inv(L) */
        cut = 0;
        while(cut < *n)
        {
            nnb = *nb;
            if(cut + nnb > *n)
            {
                nnb = *n - cut;
            }
            else
            {
                icount = 0;
                /* count negative elements, */
                i__1 = cut + nnb;
                for(i__ = cut + 1; i__ <= i__1; ++i__)
                {
                    if(ipiv[i__] < 0)
                    {
                        ++icount;
                    }
                }
                /* need a even number for a clear cut */
                if(icount % 2 == 1)
                {
                    ++nnb;
                }
            }
            /* L21 Block */
            i__1 = *n - cut - nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * work_dim1;
                    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
                }
            }
            /* L11 Block */
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = u11 + i__ + i__ * work_dim1;
                work[i__2].r = 1.f;
                work[i__2].i = 0.f; // , expr subst
                i__2 = nnb;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].r = 0.f;
                    work[i__3].i = 0.f; // , expr subst
                }
                i__2 = i__ - 1;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
                }
            }
            /* invD*L21 */
            i__ = *n - cut - nnb;
            while(i__ >= 1)
            {
                if(ipiv[cut + nnb + i__] > 0)
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        i__3 = cut + nnb + i__ + invd * work_dim1;
                        i__4 = i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[i__4].i
                                 + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                }
                else
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        u01_i_j__.r = work[i__2].r;
                        u01_i_j__.i = work[i__2].i; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        u01_ip1_j__.r = work[i__2].r;
                        u01_ip1_j__.i = work[i__2].i; // , expr subst
                        i__2 = i__ + j * work_dim1;
                        i__3 = cut + nnb + i__ + invd * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i
                                 + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i
                                 + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i
                                 + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i
                                 + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                    --i__;
                }
                --i__;
            }
            /* invD1*L11 */
            i__ = nnb;
            while(i__ >= 1)
            {
                if(ipiv[cut + i__] > 0)
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[i__4].i
                                 + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                }
                else
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        u11_i_j__.r = work[i__2].r;
                        u11_i_j__.i = work[i__2].i; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        u11_ip1_j__.r = work[i__2].r;
                        u11_ip1_j__.i = work[i__2].i; // , expr subst
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__2.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__2.i = work[i__3].r * work[i__4].i
                                 + work[i__3].i * work[i__4].r; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        q__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i * u11_ip1_j__.i;
                        q__3.i = work[i__5].r * u11_ip1_j__.i
                                 + work[i__5].i * u11_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
                        q__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * u11_i_j__.i;
                        q__2.i = work[i__3].r * u11_i_j__.i
                                 + work[i__3].i * u11_i_j__.r; // , expr subst
                        i__4 = cut + i__ - 1 + invd * work_dim1;
                        q__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i * u11_ip1_j__.i;
                        q__3.i = work[i__4].r * u11_ip1_j__.i
                                 + work[i__4].i * u11_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                    --i__;
                }
                --i__;
            }
            /* L11**H * invD1 * L11 -> L11 */
            i__1 = *n + *nb + 1;
            ctrmm_("L", uplo, "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1], lda,
                   &work[u11 + 1 + work_dim1], &i__1);
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = u11 + i__ + j * work_dim1;
                    a[i__3].r = work[i__4].r;
                    a[i__3].i = work[i__4].i; // , expr subst
                }
            }
            if(cut + nnb < *n)
            {
                /* L21**H * invD2*L21 -> A( CUT+I, CUT+J ) */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                i__3 = *n + *nb + 1;
                cgemm_("C", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + (cut + 1) * a_dim1],
                       lda, &work[work_offset], &i__2, &c_b2, &work[u11 + 1 + work_dim1], &i__3);
                /* L11 = L11**H * invD1 * L11 + U01**H * invD * U01 */
                i__1 = nnb;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = cut + i__ + (cut + j) * a_dim1;
                        i__4 = cut + i__ + (cut + j) * a_dim1;
                        i__5 = u11 + i__ + j * work_dim1;
                        q__1.r = a[i__4].r + work[i__5].r;
                        q__1.i = a[i__4].i + work[i__5].i; // , expr subst
                        a[i__3].r = q__1.r;
                        a[i__3].i = q__1.i; // , expr subst
                    }
                }
                /* L01 = L22**H * invD2 * L21 */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                ctrmm_("L", uplo, "C", "U", &i__1, &nnb, &c_b1,
                       &a[cut + nnb + 1 + (cut + nnb + 1) * a_dim1], lda, &work[work_offset],
                       &i__2);
                /* Update L21 */
                i__1 = *n - cut - nnb;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = nnb;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
                        i__4 = i__ + j * work_dim1;
                        a[i__3].r = work[i__4].r;
                        a[i__3].i = work[i__4].i; // , expr subst
                    }
                }
            }
            else
            {
                /* L11 = L11**H * invD1 * L11 */
                i__1 = nnb;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = cut + i__ + (cut + j) * a_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        a[i__3].r = work[i__4].r;
                        a[i__3].i = work[i__4].i; // , expr subst
                    }
                }
            }
            /* Next Block */
            cut += nnb;
        }
        /* Apply PERMUTATIONS P and P**T: */
        /* P * inv(L**H) * inv(D) * inv(L) * P**T. */
        /* Interchange rows and columns I and IPIV(I) in reverse order */
        /* from the formation order of IPIV vector for Lower case. */
        /* ( We can use a loop over IPIV with increment -1, */
        /* since the ABS value of IPIV(I) represents the row (column) */
        /* index of the interchange with row (column) i in both 1x1 */
        /* and 2x2 pivot cases, i.e. we don't need separate code branches */
        /* for 1x1 and 2x2 pivot cases ) */
        for(i__ = *n; i__ >= 1; --i__)
        {
            ip = (i__1 = ipiv[i__], f2c_abs(i__1));
            if(ip != i__)
            {
                if(i__ < ip)
                {
                    cheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHETRI_3X */
}
/* chetri_3x__ */
