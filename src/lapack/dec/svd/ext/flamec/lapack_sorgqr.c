/* ../netlib/sorgqr.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
/* > \brief \b SORGQR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORGQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgqr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgqr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgqr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
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
/* > SORGQR generates an M-by-N real matrix Q with orthonormal columns, */
/* > which is defined as the first N columns of a product of K elementary */
/* > reflectors of order M */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by SGEQRF. */
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
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the i-th column must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by SGEQRF in the first k columns of its array */
/* > argument A. */
/* > On exit, the M-by-N matrix Q. */
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
/* > TAU is REAL array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by SGEQRF. */
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
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
/* > For optimum performance LWORK >= N*NB, where NB is the */
/* > optimal blocksize. */
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
/* > < 0: if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int lapack_sorgqr(integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j, l, nbmin, iinfo, ib;
    extern /* Subroutine */
    int lapack_sorg2r(integer *, integer *, integer *, real *, integer *, real *, real *, integer *);
    integer nb, ki, kk, nx;
    extern /* Subroutine */
    void slarfb_(char *, char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *), xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    void slarft_(char *, char *, integer *, integer *, real *, integer *, real *, real *, integer *);
    integer ldwork, lwkopt;
    logical lquery;
    integer iws;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
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
    nb = ilaenv_(&c__1, "SORGQR", " ", m, n, k, &c_n1);
    lwkopt = fla_max(1,*n) * nb;
    work[1] = sroundup_lwork(&lwkopt);
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *n)
    {
        *info = -3;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -5;
    }
    else if (*lwork < fla_max(1,*n) && ! lquery)
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORGQR", &i__1, (ftnlen)6);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        work[1] = 1.f;
        return 0;
    }
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k)
    {
        /* Determine when to cross over from blocked to unblocked code. */
        /* Computing MAX */
        i__1 = 0;
        i__2 = ilaenv_(&c__3, "SORGQR", " ", m, n, k, &c_n1); // , expr subst
        nx = fla_max(i__1,i__2);
        if (nx < *k)
        {
            /* Determine if workspace is large enough for blocked code. */
            ldwork = *n;
            iws = ldwork * nb;
            if (*lwork < iws)
            {
                /* Not enough workspace to use optimal NB: reduce NB and */
                /* determine the minimum value of NB. */
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = ilaenv_(&c__2, "SORGQR", " ", m, n, k, &c_n1); // , expr subst
                nbmin = fla_max(i__1,i__2);
            }
        }
    }
    if (nb >= nbmin && nb < *k && nx < *k)
    {
        /* Use blocked code after the last block. */
        /* The first kk columns are handled by the block method. */
        ki = (*k - nx - 1) / nb * nb;
        /* Computing MIN */
        i__1 = *k;
        i__2 = ki + nb; // , expr subst
        kk = fla_min(i__1,i__2);
        /* Set A(1:kk,kk+1:n) to zero. */
        i__1 = *n;
        for (j = kk + 1;
                j <= i__1;
                ++j)
        {
            i__2 = kk;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] = 0.f;
                /* L10: */
            }
            /* L20: */
        }
    }
    else
    {
        kk = 0;
    }
    /* Use unblocked code for the last or only block. */
    if (kk < *n)
    {
        i__1 = *m - kk;
        i__2 = *n - kk;
        i__3 = *k - kk;
        lapack_sorg2r(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, & tau[kk + 1], &work[1], &iinfo);
    }
    if (kk > 0)
    {
        /* Use blocked code */
        i__1 = -nb;
        for (i__ = ki + 1;
                i__1 < 0 ? i__ >= 1 : i__ <= 1;
                i__ += i__1)
        {
            /* Computing MIN */
            i__2 = nb;
            i__3 = *k - i__ + 1; // , expr subst
            ib = fla_min(i__2,i__3);
            if (i__ + ib <= *n)
            {
                /* Form the triangular factor of the block reflector */
                /* H = H(i) H(i+1) . . . H(i+ib-1) */
                i__2 = *m - i__ + 1;
                slarft_("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &ldwork);
                /* Apply H to A(i:m,i+ib:n) from the left */
                i__2 = *m - i__ + 1;
                i__3 = *n - i__ - ib + 1;
                slarfb_("Left", "No transpose", "Forward", "Columnwise", & i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[ 1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, & work[ib + 1], &ldwork);
            }
            /* Apply H to rows i:m of current block */
            i__2 = *m - i__ + 1;
            lapack_sorg2r(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], & work[1], &iinfo);
            /* Set rows 1:i-1 of current block to zero */
            i__2 = i__ + ib - 1;
            for (j = i__;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ - 1;
                for (l = 1;
                        l <= i__3;
                        ++l)
                {
                    a[l + j * a_dim1] = 0.f;
                    /* L30: */
                }
                /* L40: */
            }
            /* L50: */
        }
    }
    work[1] = sroundup_lwork(&iws);
    return 0;
    /* End of SORGQR */
}
/* lapack_sorgqr */
