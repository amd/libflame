/* ../netlib/chbtrd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {{0.f}, {0.f}};
static scomplex c_b2 = {{1.f}, {0.f}};
static aocl_int64_t c__1 = 1;
/* > \brief \b CHBTRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHBTRD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbtrd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbtrd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbtrd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, */
/* WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO, VECT */
/* INTEGER INFO, KD, LDAB, LDQ, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ) */
/* COMPLEX AB( LDAB, * ), Q( LDQ, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBTRD reduces a scomplex Hermitian band matrix A to real symmetric */
/* > tridiagonal form T by a unitary similarity transformation: */
/* > Q**H * A * Q = T. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] VECT */
/* > \verbatim */
/* > VECT is CHARACTER*1 */
/* > = 'N': do not form Q;
 */
/* > = 'V': form Q;
 */
/* > = 'U': update a matrix X, by forming X*Q. */
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
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > On entry, the upper or lower triangle of the Hermitian band */
/* > matrix A, stored in the first KD+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+kd). */
/* > On exit, the diagonal elements of AB are overwritten by the */
/* > diagonal elements of the tridiagonal matrix T;
if KD > 0, the */
/* > elements on the first superdiagonal (if UPLO = 'U') or the */
/* > first subdiagonal (if UPLO = 'L') are overwritten by the */
/* > off-diagonal elements of T;
the rest of AB is overwritten by */
/* > values generated during the reduction. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix T: */
/* > E(i) = T(i,i+1) if UPLO = 'U';
E(i) = T(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ,N) */
/* > On entry, if VECT = 'U', then Q must contain an N-by-N */
/* > matrix X;
if VECT = 'N' or 'V', then Q need not be set. */
/* > */
/* > On exit: */
/* > if VECT = 'V', Q contains the N-by-N unitary matrix Q;
 */
/* > if VECT = 'U', Q contains the product X*Q;
 */
/* > if VECT = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. */
/* > LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
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
/* > Modified by Linda Kaufman, Bell Labs. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chbtrd_(char *vect, char *uplo, aocl_int_t *n, aocl_int_t *kd, scomplex *ab, aocl_int_t *ldab,
             real *d__, real *e, scomplex *q, aocl_int_t *ldq, scomplex *work, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chbtrd(vect, uplo, n, kd, ab, ldab, d__, e, q, ldq, work, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t kd_64 = *kd;
    aocl_int64_t ldab_64 = *ldab;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chbtrd(vect, uplo, &n_64, &kd_64, ab, &ldab_64, d__, e, q, &ldq_64, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_chbtrd(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab,
                        aocl_int64_t *ldab, real *d__, real *e, scomplex *q, aocl_int64_t *ldq,
                        scomplex *work, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "chbtrd inputs: vect %c, uplo %c, n %lld, kd %lld, ldab %lld, ldq %lld",
             *vect, *uplo, *n, *kd, *ldab, *ldq);
#else
    snprintf(buffer, 256, "chbtrd inputs: vect %c, uplo %c, n %d, kd %d, ldab %d, ldq %d", *vect,
             *uplo, *n, *kd, *ldab, *ldq);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t ab_dim1, ab_offset, q_dim1, q_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1;
    scomplex q__1;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    double c_abs(scomplex *);
    /* Local variables */
    aocl_int64_t i__, j, k, l;
    scomplex t;
    aocl_int64_t i2, j1, j2, nq, nr, kd1, ibl, iqb, kdn, jin, nrt, kdm1, inca, jend, lend, jinc;
    real abst;
    aocl_int64_t incx, last;
    scomplex temp;
    aocl_int64_t j1end, j1inc;
    aocl_int64_t iqend;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical initq, wantq, upper;
    aocl_int64_t iqaend;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --d__;
    --e;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;
    /* Function Body */
    initq = lsame_(vect, "V", 1, 1);
    wantq = initq || lsame_(vect, "U", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    kd1 = *kd + 1;
    kdm1 = *kd - 1;
    incx = *ldab - 1;
    iqend = 1;
    *info = 0;
    if(!wantq && !lsame_(vect, "N", 1, 1))
    {
        *info = -1;
    }
    else if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*kd < 0)
    {
        *info = -4;
    }
    else if(*ldab < kd1)
    {
        *info = -6;
    }
    else if(*ldq < fla_max(1, *n) && wantq)
    {
        *info = -10;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CHBTRD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Initialize Q to the unit matrix, if needed */
    if(initq)
    {
        aocl_lapack_claset("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq);
    }
    /* Wherever possible, plane rotations are generated and applied in */
    /* vector operations of length NR over the index set J1:J2:KD1. */
    /* The real cosines and scomplex sines of the plane rotations are */
    /* stored in the arrays D and WORK. */
    inca = kd1 * *ldab;
    /* Computing MIN */
    i__1 = *n - 1;
    kdn = fla_min(i__1, *kd);
    if(upper)
    {
        if(*kd > 1)
        {
            /* Reduce to scomplex Hermitian tridiagonal form, working with */
            /* the upper triangle */
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;
            i__1 = kd1 + ab_dim1;
            i__2 = kd1 + ab_dim1;
            r__1 = ab[i__2].r;
            ab[i__1].r = r__1;
            ab[i__1].i = 0.f; // , expr subst
            i__1 = *n - 2;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                /* Reduce i-th row of matrix to tridiagonal form */
                for(k = kdn + 1; k >= 2; --k)
                {
                    j1 += kdn;
                    j2 += kdn;
                    if(nr > 0)
                    {
                        /* generate plane rotations to annihilate nonzero */
                        /* elements which have been created outside the band */
                        aocl_lapack_clargv(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &inca, &work[j1], &kd1,
                                           &d__[j1], &kd1);
                        /* apply rotations from the right */
                        /* Dependent on the the number of diagonals either */
                        /* CLARTV or CROT is used */
                        if(nr >= (*kd << 1) - 1)
                        {
                            i__2 = *kd - 1;
                            for(l = 1; l <= i__2; ++l)
                            {
                                aocl_lapack_clartv(&nr, &ab[l + 1 + (j1 - 1) * ab_dim1], &inca,
                                                   &ab[l + j1 * ab_dim1], &inca, &d__[j1],
                                                   &work[j1], &kd1);
                                /* L10: */
                            }
                        }
                        else
                        {
                            jend = j1 + (nr - 1) * kd1;
                            i__2 = jend;
                            i__3 = kd1;
                            for(jinc = j1; i__3 < 0 ? jinc >= i__2 : jinc <= i__2; jinc += i__3)
                            {
                                aocl_lapack_crot(&kdm1, &ab[(jinc - 1) * ab_dim1 + 2], &c__1,
                                                 &ab[jinc * ab_dim1 + 1], &c__1, &d__[jinc],
                                                 &work[jinc]);
                                /* L20: */
                            }
                        }
                    }
                    if(k > 2)
                    {
                        if(k <= *n - i__ + 1)
                        {
                            /* generate plane rotation to annihilate a(i,i+k-1) */
                            /* within the band */
                            clartg_(&ab[*kd - k + 3 + (i__ + k - 2) * ab_dim1],
                                    &ab[*kd - k + 2 + (i__ + k - 1) * ab_dim1], &d__[i__ + k - 1],
                                    &work[i__ + k - 1], &temp);
                            i__3 = *kd - k + 3 + (i__ + k - 2) * ab_dim1;
                            ab[i__3].r = temp.r;
                            ab[i__3].i = temp.i; // , expr subst
                            /* apply rotation from the right */
                            i__3 = k - 3;
                            aocl_lapack_crot(&i__3, &ab[*kd - k + 4 + (i__ + k - 2) * ab_dim1],
                                             &c__1, &ab[*kd - k + 3 + (i__ + k - 1) * ab_dim1],
                                             &c__1, &d__[i__ + k - 1], &work[i__ + k - 1]);
                        }
                        ++nr;
                        j1 = j1 - kdn - 1;
                    }
                    /* apply plane rotations from both sides to diagonal */
                    /* blocks */
                    if(nr > 0)
                    {
                        aocl_lapack_clar2v(&nr, &ab[kd1 + (j1 - 1) * ab_dim1],
                                           &ab[kd1 + j1 * ab_dim1], &ab[*kd + j1 * ab_dim1], &inca,
                                           &d__[j1], &work[j1], &kd1);
                    }
                    /* apply plane rotations from the left */
                    if(nr > 0)
                    {
                        aocl_lapack_clacgv(&nr, &work[j1], &kd1);
                        if((*kd << 1) - 1 < nr)
                        {
                            /* Dependent on the the number of diagonals either */
                            /* CLARTV or CROT is used */
                            i__3 = *kd - 1;
                            for(l = 1; l <= i__3; ++l)
                            {
                                if(j2 + l > *n)
                                {
                                    nrt = nr - 1;
                                }
                                else
                                {
                                    nrt = nr;
                                }
                                if(nrt > 0)
                                {
                                    aocl_lapack_clartv(&nrt, &ab[*kd - l + (j1 + l) * ab_dim1],
                                                       &inca, &ab[*kd - l + 1 + (j1 + l) * ab_dim1],
                                                       &inca, &d__[j1], &work[j1], &kd1);
                                }
                                /* L30: */
                            }
                        }
                        else
                        {
                            j1end = j1 + kd1 * (nr - 2);
                            if(j1end >= j1)
                            {
                                i__3 = j1end;
                                i__2 = kd1;
                                for(jin = j1; i__2 < 0 ? jin >= i__3 : jin <= i__3; jin += i__2)
                                {
                                    i__4 = *kd - 1;
                                    aocl_lapack_crot(&i__4, &ab[*kd - 1 + (jin + 1) * ab_dim1],
                                                     &incx, &ab[*kd + (jin + 1) * ab_dim1], &incx,
                                                     &d__[jin], &work[jin]);
                                    /* L40: */
                                }
                            }
                            /* Computing MIN */
                            i__2 = kdm1;
                            i__3 = *n - j2; // , expr subst
                            lend = fla_min(i__2, i__3);
                            last = j1end + kd1;
                            if(lend > 0)
                            {
                                aocl_lapack_crot(&lend, &ab[*kd - 1 + (last + 1) * ab_dim1], &incx,
                                                 &ab[*kd + (last + 1) * ab_dim1], &incx, &d__[last],
                                                 &work[last]);
                            }
                        }
                    }
                    if(wantq)
                    {
                        /* accumulate product of plane rotations in Q */
                        if(initq)
                        {
                            /* take advantage of the fact that Q was */
                            /* initially the Identity matrix */
                            iqend = fla_max(iqend, j2);
                            /* Computing MAX */
                            i__2 = 0;
                            i__3 = k - 3; // , expr subst
                            i2 = fla_max(i__2, i__3);
                            iqaend = i__ * *kd + 1;
                            if(k == 2)
                            {
                                iqaend += *kd;
                            }
                            iqaend = fla_min(iqaend, iqend);
                            i__2 = j2;
                            i__3 = kd1;
                            for(j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3)
                            {
                                ibl = i__ - i2 / kdm1;
                                ++i2;
                                /* Computing MAX */
                                i__4 = 1;
                                i__5 = j - ibl; // , expr subst
                                iqb = fla_max(i__4, i__5);
                                nq = iqaend + 1 - iqb;
                                /* Computing MIN */
                                i__4 = iqaend + *kd;
                                iqaend = fla_min(i__4, iqend);
                                r_cnjg(&q__1, &work[j]);
                                aocl_lapack_crot(&nq, &q[iqb + (j - 1) * q_dim1], &c__1,
                                                 &q[iqb + j * q_dim1], &c__1, &d__[j], &q__1);
                                /* L50: */
                            }
                        }
                        else
                        {
                            i__3 = j2;
                            i__2 = kd1;
                            for(j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2)
                            {
                                r_cnjg(&q__1, &work[j]);
                                aocl_lapack_crot(n, &q[(j - 1) * q_dim1 + 1], &c__1,
                                                 &q[j * q_dim1 + 1], &c__1, &d__[j], &q__1);
                                /* L60: */
                            }
                        }
                    }
                    if(j2 + kdn > *n)
                    {
                        /* adjust J2 to keep within the bounds of the matrix */
                        --nr;
                        j2 = j2 - kdn - 1;
                    }
                    i__2 = j2;
                    i__3 = kd1;
                    for(j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3)
                    {
                        /* create nonzero element a(j-1,j+kd) outside the band */
                        /* and store it in WORK */
                        i__4 = j + *kd;
                        i__5 = j;
                        i__6 = (j + *kd) * ab_dim1 + 1;
                        q__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6].i;
                        q__1.i
                            = work[i__5].r * ab[i__6].i + work[i__5].i * ab[i__6].r; // , expr subst
                        work[i__4].r = q__1.r;
                        work[i__4].i = q__1.i; // , expr subst
                        i__4 = (j + *kd) * ab_dim1 + 1;
                        i__5 = j;
                        i__6 = (j + *kd) * ab_dim1 + 1;
                        q__1.r = d__[i__5] * ab[i__6].r;
                        q__1.i = d__[i__5] * ab[i__6].i; // , expr subst
                        ab[i__4].r = q__1.r;
                        ab[i__4].i = q__1.i; // , expr subst
                        /* L70: */
                    }
                    /* L80: */
                }
                /* L90: */
            }
        }
        if(*kd > 0)
        {
            /* make off-diagonal elements real and copy them to E */
            i__1 = *n - 1;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__3 = *kd + (i__ + 1) * ab_dim1;
                t.r = ab[i__3].r;
                t.i = ab[i__3].i; // , expr subst
                abst = c_abs(&t);
                i__3 = *kd + (i__ + 1) * ab_dim1;
                ab[i__3].r = abst;
                ab[i__3].i = 0.f; // , expr subst
                e[i__] = abst;
                if(abst != 0.f)
                {
                    q__1.r = t.r / abst;
                    q__1.i = t.i / abst; // , expr subst
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                }
                else
                {
                    t.r = 1.f;
                    t.i = 0.f; // , expr subst
                }
                if(i__ < *n - 1)
                {
                    i__3 = *kd + (i__ + 2) * ab_dim1;
                    i__2 = *kd + (i__ + 2) * ab_dim1;
                    q__1.r = ab[i__2].r * t.r - ab[i__2].i * t.i;
                    q__1.i = ab[i__2].r * t.i + ab[i__2].i * t.r; // , expr subst
                    ab[i__3].r = q__1.r;
                    ab[i__3].i = q__1.i; // , expr subst
                }
                if(wantq)
                {
                    r_cnjg(&q__1, &t);
                    aocl_blas_cscal(n, &q__1, &q[(i__ + 1) * q_dim1 + 1], &c__1);
                }
                /* L100: */
            }
        }
        else
        {
            /* set E to zero if original matrix was diagonal */
            i__1 = *n - 1;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                e[i__] = 0.f;
                /* L110: */
            }
        }
        /* copy diagonal elements to D */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__3 = i__;
            i__2 = kd1 + i__ * ab_dim1;
            d__[i__3] = ab[i__2].r;
            /* L120: */
        }
    }
    else
    {
        if(*kd > 1)
        {
            /* Reduce to scomplex Hermitian tridiagonal form, working with */
            /* the lower triangle */
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;
            i__1 = ab_dim1 + 1;
            i__3 = ab_dim1 + 1;
            r__1 = ab[i__3].r;
            ab[i__1].r = r__1;
            ab[i__1].i = 0.f; // , expr subst
            i__1 = *n - 2;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                /* Reduce i-th column of matrix to tridiagonal form */
                for(k = kdn + 1; k >= 2; --k)
                {
                    j1 += kdn;
                    j2 += kdn;
                    if(nr > 0)
                    {
                        /* generate plane rotations to annihilate nonzero */
                        /* elements which have been created outside the band */
                        aocl_lapack_clargv(&nr, &ab[kd1 + (j1 - kd1) * ab_dim1], &inca, &work[j1],
                                           &kd1, &d__[j1], &kd1);
                        /* apply plane rotations from one side */
                        /* Dependent on the the number of diagonals either */
                        /* CLARTV or CROT is used */
                        if(nr > (*kd << 1) - 1)
                        {
                            i__3 = *kd - 1;
                            for(l = 1; l <= i__3; ++l)
                            {
                                aocl_lapack_clartv(&nr, &ab[kd1 - l + (j1 - kd1 + l) * ab_dim1],
                                                   &inca,
                                                   &ab[kd1 - l + 1 + (j1 - kd1 + l) * ab_dim1],
                                                   &inca, &d__[j1], &work[j1], &kd1);
                                /* L130: */
                            }
                        }
                        else
                        {
                            jend = j1 + kd1 * (nr - 1);
                            i__3 = jend;
                            i__2 = kd1;
                            for(jinc = j1; i__2 < 0 ? jinc >= i__3 : jinc <= i__3; jinc += i__2)
                            {
                                aocl_lapack_crot(&kdm1, &ab[*kd + (jinc - *kd) * ab_dim1], &incx,
                                                 &ab[kd1 + (jinc - *kd) * ab_dim1], &incx,
                                                 &d__[jinc], &work[jinc]);
                                /* L140: */
                            }
                        }
                    }
                    if(k > 2)
                    {
                        if(k <= *n - i__ + 1)
                        {
                            /* generate plane rotation to annihilate a(i+k-1,i) */
                            /* within the band */
                            clartg_(&ab[k - 1 + i__ * ab_dim1], &ab[k + i__ * ab_dim1],
                                    &d__[i__ + k - 1], &work[i__ + k - 1], &temp);
                            i__2 = k - 1 + i__ * ab_dim1;
                            ab[i__2].r = temp.r;
                            ab[i__2].i = temp.i; // , expr subst
                            /* apply rotation from the left */
                            i__2 = k - 3;
                            i__3 = *ldab - 1;
                            i__4 = *ldab - 1;
                            aocl_lapack_crot(&i__2, &ab[k - 2 + (i__ + 1) * ab_dim1], &i__3,
                                             &ab[k - 1 + (i__ + 1) * ab_dim1], &i__4,
                                             &d__[i__ + k - 1], &work[i__ + k - 1]);
                        }
                        ++nr;
                        j1 = j1 - kdn - 1;
                    }
                    /* apply plane rotations from both sides to diagonal */
                    /* blocks */
                    if(nr > 0)
                    {
                        aocl_lapack_clar2v(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &ab[j1 * ab_dim1 + 1],
                                           &ab[(j1 - 1) * ab_dim1 + 2], &inca, &d__[j1], &work[j1],
                                           &kd1);
                    }
                    /* apply plane rotations from the right */
                    /* Dependent on the the number of diagonals either */
                    /* CLARTV or CROT is used */
                    if(nr > 0)
                    {
                        aocl_lapack_clacgv(&nr, &work[j1], &kd1);
                        if(nr > (*kd << 1) - 1)
                        {
                            i__2 = *kd - 1;
                            for(l = 1; l <= i__2; ++l)
                            {
                                if(j2 + l > *n)
                                {
                                    nrt = nr - 1;
                                }
                                else
                                {
                                    nrt = nr;
                                }
                                if(nrt > 0)
                                {
                                    aocl_lapack_clartv(&nrt, &ab[l + 2 + (j1 - 1) * ab_dim1], &inca,
                                                       &ab[l + 1 + j1 * ab_dim1], &inca, &d__[j1],
                                                       &work[j1], &kd1);
                                }
                                /* L150: */
                            }
                        }
                        else
                        {
                            j1end = j1 + kd1 * (nr - 2);
                            if(j1end >= j1)
                            {
                                i__2 = j1end;
                                i__3 = kd1;
                                for(j1inc = j1; i__3 < 0 ? j1inc >= i__2 : j1inc <= i__2;
                                    j1inc += i__3)
                                {
                                    aocl_lapack_crot(&kdm1, &ab[(j1inc - 1) * ab_dim1 + 3], &c__1,
                                                     &ab[j1inc * ab_dim1 + 2], &c__1, &d__[j1inc],
                                                     &work[j1inc]);
                                    /* L160: */
                                }
                            }
                            /* Computing MIN */
                            i__3 = kdm1;
                            i__2 = *n - j2; // , expr subst
                            lend = fla_min(i__3, i__2);
                            last = j1end + kd1;
                            if(lend > 0)
                            {
                                aocl_lapack_crot(&lend, &ab[(last - 1) * ab_dim1 + 3], &c__1,
                                                 &ab[last * ab_dim1 + 2], &c__1, &d__[last],
                                                 &work[last]);
                            }
                        }
                    }
                    if(wantq)
                    {
                        /* accumulate product of plane rotations in Q */
                        if(initq)
                        {
                            /* take advantage of the fact that Q was */
                            /* initially the Identity matrix */
                            iqend = fla_max(iqend, j2);
                            /* Computing MAX */
                            i__3 = 0;
                            i__2 = k - 3; // , expr subst
                            i2 = fla_max(i__3, i__2);
                            iqaend = i__ * *kd + 1;
                            if(k == 2)
                            {
                                iqaend += *kd;
                            }
                            iqaend = fla_min(iqaend, iqend);
                            i__3 = j2;
                            i__2 = kd1;
                            for(j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2)
                            {
                                ibl = i__ - i2 / kdm1;
                                ++i2;
                                /* Computing MAX */
                                i__4 = 1;
                                i__5 = j - ibl; // , expr subst
                                iqb = fla_max(i__4, i__5);
                                nq = iqaend + 1 - iqb;
                                /* Computing MIN */
                                i__4 = iqaend + *kd;
                                iqaend = fla_min(i__4, iqend);
                                aocl_lapack_crot(&nq, &q[iqb + (j - 1) * q_dim1], &c__1,
                                                 &q[iqb + j * q_dim1], &c__1, &d__[j], &work[j]);
                                /* L170: */
                            }
                        }
                        else
                        {
                            i__2 = j2;
                            i__3 = kd1;
                            for(j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3)
                            {
                                aocl_lapack_crot(n, &q[(j - 1) * q_dim1 + 1], &c__1,
                                                 &q[j * q_dim1 + 1], &c__1, &d__[j], &work[j]);
                                /* L180: */
                            }
                        }
                    }
                    if(j2 + kdn > *n)
                    {
                        /* adjust J2 to keep within the bounds of the matrix */
                        --nr;
                        j2 = j2 - kdn - 1;
                    }
                    i__3 = j2;
                    i__2 = kd1;
                    for(j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2)
                    {
                        /* create nonzero element a(j+kd,j-1) outside the */
                        /* band and store it in WORK */
                        i__4 = j + *kd;
                        i__5 = j;
                        i__6 = kd1 + j * ab_dim1;
                        q__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6].i;
                        q__1.i
                            = work[i__5].r * ab[i__6].i + work[i__5].i * ab[i__6].r; // , expr subst
                        work[i__4].r = q__1.r;
                        work[i__4].i = q__1.i; // , expr subst
                        i__4 = kd1 + j * ab_dim1;
                        i__5 = j;
                        i__6 = kd1 + j * ab_dim1;
                        q__1.r = d__[i__5] * ab[i__6].r;
                        q__1.i = d__[i__5] * ab[i__6].i; // , expr subst
                        ab[i__4].r = q__1.r;
                        ab[i__4].i = q__1.i; // , expr subst
                        /* L190: */
                    }
                    /* L200: */
                }
                /* L210: */
            }
        }
        if(*kd > 0)
        {
            /* make off-diagonal elements real and copy them to E */
            i__1 = *n - 1;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__ * ab_dim1 + 2;
                t.r = ab[i__2].r;
                t.i = ab[i__2].i; // , expr subst
                abst = c_abs(&t);
                i__2 = i__ * ab_dim1 + 2;
                ab[i__2].r = abst;
                ab[i__2].i = 0.f; // , expr subst
                e[i__] = abst;
                if(abst != 0.f)
                {
                    q__1.r = t.r / abst;
                    q__1.i = t.i / abst; // , expr subst
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                }
                else
                {
                    t.r = 1.f;
                    t.i = 0.f; // , expr subst
                }
                if(i__ < *n - 1)
                {
                    i__2 = (i__ + 1) * ab_dim1 + 2;
                    i__3 = (i__ + 1) * ab_dim1 + 2;
                    q__1.r = ab[i__3].r * t.r - ab[i__3].i * t.i;
                    q__1.i = ab[i__3].r * t.i + ab[i__3].i * t.r; // , expr subst
                    ab[i__2].r = q__1.r;
                    ab[i__2].i = q__1.i; // , expr subst
                }
                if(wantq)
                {
                    aocl_blas_cscal(n, &t, &q[(i__ + 1) * q_dim1 + 1], &c__1);
                }
                /* L220: */
            }
        }
        else
        {
            /* set E to zero if original matrix was diagonal */
            i__1 = *n - 1;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                e[i__] = 0.f;
                /* L230: */
            }
        }
        /* copy diagonal elements to D */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = i__;
            i__3 = i__ * ab_dim1 + 1;
            d__[i__2] = ab[i__3].r;
            /* L240: */
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHBTRD */
}
/* chbtrd_ */
