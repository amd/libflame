/* ../netlib/dorg2r.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
    Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DORG2R generates all or part of the orthogonal matrix Q from a QR factorization
 * determined by s geqrf (unblocked algorithm). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORG2R + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORG2R generates an m by n real matrix Q with orthonormal columns, */
/* > which is defined as the first n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by DGEQRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix Q. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix Q. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines the */
/* > matrix Q. N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the i-th column must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by DGEQRF in the first k columns of its array */
/* > argument A. */
/* > On exit, the m-by-n matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The first dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by DGEQRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void dorg2r_fla(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau,
                doublereal *work, integer *info)
{
    extern fla_context fla_global_context;
    extern void dorg2r_fla_opt(integer * m, integer * n, integer * k, doublereal * a, integer * lda,
                               doublereal * tau, doublereal * work, integer * info);
    extern void dorg2r_fla_native(integer * m, integer * n, integer * k, doublereal * a,
                                  integer * lda, doublereal * tau, doublereal * work,
                                  integer * info);

    /* Initialize global context data */
    aocl_fla_init();

#if FLA_ENABLE_AMD_OPT
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        dorg2r_fla_opt(m, n, k, a, lda, tau, work, info);
    }
    else
    {
        dorg2r_fla_native(m, n, k, a, lda, tau, work, info);
    }
#else
    dorg2r_fla_native(m, n, k, a, lda, tau, work, info);
#endif

    return;
}

#if FLA_ENABLE_AMD_OPT
void dorg2r_fla_opt(integer *m, integer *n, integer *k, doublereal *a, integer *lda,
                    doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    /* Local variables */
    integer i__, j, l;
    extern /* Subroutine */
        void
        fla_dscal(integer *, doublereal *, doublereal *, integer *),
        dscal_(integer *, doublereal *, doublereal *, integer *),
        dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
               integer *, doublereal *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if(*k < 0 || *k > *n)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORG2R", &i__1, (ftnlen)6);
        return;
    }
    /* Quick return if possible */
    if(*n <= 0)
    {
        return;
    }
    /* Initialise columns k+1:n to columns of the unit matrix */
    i__1 = *n;
    for(j = *k + 1; j <= i__1; ++j)
    {
        i__2 = *m;
        for(l = 1; l <= i__2; ++l)
        {
            a[l + j * a_dim1] = 0.;
            /* L10: */
        }
        a[j + j * a_dim1] = 1.;
        /* L20: */
    }

    for(i__ = *k; i__ >= 1; --i__)
    {
        /* Apply H(i) to A(i:m,i:n) from the left */
        if(i__ < *n)
        {
            a[i__ + i__ * a_dim1] = 1.;
            i__1 = *m - i__ + 1;
            i__2 = *n - i__;
            dlarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[i__],
                   &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
        }

        /* Inline DSCAL for small size */
        i__1 = *m - i__;
        d__1 = -tau[i__];

        if(i__ < *m)
        {
            fla_dscal(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
        }
        else
        {
            dscal_(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
        }

        a[i__ + i__ * a_dim1] = 1. - tau[i__];
        /* Set A(1:i-1,i) to zero */
        i__1 = i__ - 1;
        for(l = 1; l <= i__1; ++l)
        {
            a[l + i__ * a_dim1] = 0.;
            /* L30: */
        }
        /* L40: */
    }
    return;
    /* End of DORG2R */
}
#endif

void dorg2r_fla_native(integer *m, integer *n, integer *k, doublereal *a, integer *lda,
                       doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    /* Local variables */
    integer i__, j, l;
    extern /* Subroutine */
        void
        dscal_(integer *, doublereal *, doublereal *, integer *),
        dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
               integer *, doublereal *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if(*k < 0 || *k > *n)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORG2R", &i__1, (ftnlen)6);
        return;
    }
    /* Quick return if possible */
    if(*n <= 0)
    {
        return;
    }
    /* Initialise columns k+1:n to columns of the unit matrix */
    i__1 = *n;
    for(j = *k + 1; j <= i__1; ++j)
    {
        i__2 = *m;
        for(l = 1; l <= i__2; ++l)
        {
            a[l + j * a_dim1] = 0.;
            /* L10: */
        }
        a[j + j * a_dim1] = 1.;
        /* L20: */
    }

    for(i__ = *k; i__ >= 1; --i__)
    {
        /* Apply H(i) to A(i:m,i:n) from the left */
        if(i__ < *n)
        {
            a[i__ + i__ * a_dim1] = 1.;
            i__1 = *m - i__ + 1;
            i__2 = *n - i__;
            dlarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[i__],
                   &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
        }

        if(i__ < *m)
        {
            i__1 = *m - i__;
            d__1 = -tau[i__];
            dscal_(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
        }
        a[i__ + i__ * a_dim1] = 1. - tau[i__];
        /* Set A(1:i-1,i) to zero */
        i__1 = i__ - 1;
        for(l = 1; l <= i__1; ++l)
        {
            a[l + i__ * a_dim1] = 0.;
            /* L30: */
        }
        /* L40: */
    }
    return;
    /* End of DORG2R */
}
/* dorg2r_ */
