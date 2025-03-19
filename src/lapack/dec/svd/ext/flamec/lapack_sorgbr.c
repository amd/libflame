/* ../netlib/sorgbr.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
/* > \brief \b SORGBR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORGBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgbr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgbr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgbr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER VECT */
/* INTEGER INFO, K, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORGBR generates one of the real orthogonal matrices Q or P**T */
/* > determined by SGEBRD when reducing a real matrix A to bidiagonal */
/* > form: A = Q * B * P**T. Q and P**T are defined as products of */
/* > elementary reflectors H(i) or G(i) respectively. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q */
/* > is of order M: */
/* > if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n */
/* > columns of Q, where m >= n >= k;
*/
/* > if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an */
/* > M-by-M matrix. */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T */
/* > is of order N: */
/* > if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m */
/* > rows of P**T, where n >= m >= k;
*/
/* > if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as */
/* > an N-by-N matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] VECT */
/* > \verbatim */
/* > VECT is CHARACTER*1 */
/* > Specifies whether the matrix Q or the matrix P**T is */
/* > required, as defined in the transformation applied by SGEBRD: */
/* > = 'Q': generate Q;
*/
/* > = 'P': generate P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix Q or P**T to be returned. */
/* > M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix Q or P**T to be returned. */
/* > N >= 0. */
/* > If VECT = 'Q', M >= N >= fla_min(M,K);
*/
/* > if VECT = 'P', N >= M >= fla_min(N,K). */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > If VECT = 'Q', the number of columns in the original M-by-K */
/* > matrix reduced by SGEBRD. */
/* > If VECT = 'P', the number of rows in the original K-by-N */
/* > matrix reduced by SGEBRD. */
/* > K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the vectors which define the elementary reflectors, */
/* > as returned by SGEBRD. */
/* > On exit, the M-by-N matrix Q or P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is REAL array, dimension */
/* > (fla_min(M,K)) if VECT = 'Q' */
/* > (fla_min(N,K)) if VECT = 'P' */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i) or G(i), which determines Q or P**T, as */
/* > returned by SGEBRD in its array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,fla_min(M,N)). */
/* > For optimum performance LWORK >= fla_min(M,N)*NB, where NB */
/* > is the optimal blocksize. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
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
/* > \ingroup realGBcomputational */
/* ===================================================================== */
/* Subroutine */
int lapack_sorgbr(char *vect, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    logical wantq;
    integer mn;
    extern /* Subroutine */
    void xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern int lapack_sorglq( integer *, integer *, integer *, real *, integer *, real *, real *, integer *, integer *), lapack_sorgqr(integer *, integer *, integer *, real *, integer *, real *, real *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    extern real sroundup_lwork(integer *);
    /* -- LAPACK computational routine -- */
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
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    wantq = lsame_(vect, "Q", 1, 1);
    mn = fla_min(*m,*n);
    lquery = *lwork == -1;
    if (! wantq && ! lsame_(vect, "P", 1, 1))
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*n < 0 || wantq && (*n > *m || *n < fla_min(*m,*k)) || ! wantq && ( *m > *n || *m < fla_min(*n,*k)))
    {
        *info = -3;
    }
    else if (*k < 0)
    {
        *info = -4;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -6;
    }
    else if (*lwork < fla_max(1,mn) && ! lquery)
    {
        *info = -9;
    }
    if (*info == 0)
    {
        work[1] = 1.f;
        if (wantq)
        {
            if (*m >= *k)
            {
                lapack_sorgqr(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*m > 1)
                {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    lapack_sorgqr(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], & work[1], &c_n1, &iinfo);
                }
            }
        }
        else
        {
            if (*k < *n)
            {
                lapack_sorglq(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if (*n > 1)
                {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    lapack_sorglq(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], & work[1], &c_n1, &iinfo);
                }
            }
        }
        lwkopt = (integer)work[1];
        lwkopt = fla_max(lwkopt,mn);
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORGBR", &i__1, (ftnlen)6);
        return 0;
    }
    else if (lquery)
    {
        work[1] = sroundup_lwork(&lwkopt);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        work[1] = 1.f;
        return 0;
    }
    if (wantq)
    {
        /* Form Q, determined by a call to SGEBRD to reduce an m-by-k */
        /* matrix */
        if (*m >= *k)
        {
            /* If m >= k, assume m >= n >= k */
            lapack_sorgqr(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, & iinfo);
        }
        else
        {
            /* If m < k, assume m = n */
            /* Shift the vectors which define the elementary reflectors one */
            /* column to the right, and set the first row and column of Q */
            /* to those of the unit matrix */
            for (j = *m;
                    j >= 2;
                    --j)
            {
                a[j * a_dim1 + 1] = 0.f;
                i__1 = *m;
                for (i__ = j + 1;
                        i__ <= i__1;
                        ++i__)
                {
                    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
                    /* L10: */
                }
                /* L20: */
            }
            a[a_dim1 + 1] = 1.f;
            i__1 = *m;
            for (i__ = 2;
                    i__ <= i__1;
                    ++i__)
            {
                a[i__ + a_dim1] = 0.f;
                /* L30: */
            }
            if (*m > 1)
            {
                /* Form Q(2:m,2:m) */
                i__1 = *m - 1;
                i__2 = *m - 1;
                i__3 = *m - 1;
                lapack_sorgqr(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[ 1], &work[1], lwork, &iinfo);
            }
        }
    }
    else
    {
        /* Form P**T, determined by a call to SGEBRD to reduce a k-by-n */
        /* matrix */
        if (*k < *n)
        {
            /* If k < n, assume k <= m <= n */
            lapack_sorglq(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, & iinfo);
        }
        else
        {
            /* If k >= n, assume m = n */
            /* Shift the vectors which define the elementary reflectors one */
            /* row downward, and set the first row and column of P**T to */
            /* those of the unit matrix */
            a[a_dim1 + 1] = 1.f;
            i__1 = *n;
            for (i__ = 2;
                    i__ <= i__1;
                    ++i__)
            {
                a[i__ + a_dim1] = 0.f;
                /* L40: */
            }
            i__1 = *n;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                for (i__ = j - 1;
                        i__ >= 2;
                        --i__)
                {
                    a[i__ + j * a_dim1] = a[i__ - 1 + j * a_dim1];
                    /* L50: */
                }
                a[j * a_dim1 + 1] = 0.f;
                /* L60: */
            }
            if (*n > 1)
            {
                /* Form P**T(2:n,2:n) */
                i__1 = *n - 1;
                i__2 = *n - 1;
                i__3 = *n - 1;
                lapack_sorglq(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[ 1], &work[1], lwork, &iinfo);
            }
        }
    }
    work[1] = sroundup_lwork(&lwkopt);
    return 0;
    /* End of SORGBR */
}
/* lapack_sorgbr */
