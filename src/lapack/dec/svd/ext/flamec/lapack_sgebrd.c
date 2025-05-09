/* ../netlib/sgebrd.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
/* > \brief \b SGEBRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGEBRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgebrd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgebrd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgebrd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), D( * ), E( * ), TAUP( * ), */
/* $ TAUQ( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEBRD reduces a general real M-by-N matrix A to upper or lower */
/* > bidiagonal form B by an orthogonal transformation: Q**T * A * P = B. */
/* > */
/* > If m >= n, B is upper bidiagonal;
if m < n, B is lower bidiagonal. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows in the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns in the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N general matrix to be reduced. */
/* > On exit, */
/* > if m >= n, the diagonal and the first superdiagonal are */
/* > overwritten with the upper bidiagonal matrix B;
the */
/* > elements below the diagonal, with the array TAUQ, represent */
/* > the orthogonal matrix Q as a product of elementary */
/* > reflectors, and the elements above the first superdiagonal, */
/* > with the array TAUP, represent the orthogonal matrix P as */
/* > a product of elementary reflectors;
*/
/* > if m < n, the diagonal and the first subdiagonal are */
/* > overwritten with the lower bidiagonal matrix B;
the */
/* > elements below the first subdiagonal, with the array TAUQ, */
/* > represent the orthogonal matrix Q as a product of */
/* > elementary reflectors, and the elements above the diagonal, */
/* > with the array TAUP, represent the orthogonal matrix P as */
/* > a product of elementary reflectors. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (fla_min(M,N)) */
/* > The diagonal elements of the bidiagonal matrix B: */
/* > D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (fla_min(M,N)-1) */
/* > The off-diagonal elements of the bidiagonal matrix B: */
/* > if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*/
/* > if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* > TAUQ is REAL array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the orthogonal matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* > TAUP is REAL array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the orthogonal matrix P. See Further Details. */
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
/* > The length of the array WORK. LWORK >= fla_max(1,M,N). */
/* > For optimum performance LWORK >= (M+N)*NB, where NB */
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
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup realGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrices Q and P are represented as products of elementary */
/* > reflectors: */
/* > */
/* > If m >= n, */
/* > */
/* > Q = H(1) H(2) . . . H(n) and P = G(1) G(2) . . . G(n-1) */
/* > */
/* > Each H(i) and G(i) has the form: */
/* > */
/* > H(i) = I - tauq * v * v**T and G(i) = I - taup * u * u**T */
/* > */
/* > where tauq and taup are real scalars, and v and u are real vectors;
*/
/* > v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
*/
/* > u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
*/
/* > tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* > If m < n, */
/* > */
/* > Q = H(1) H(2) . . . H(m-1) and P = G(1) G(2) . . . G(m) */
/* > */
/* > Each H(i) and G(i) has the form: */
/* > */
/* > H(i) = I - tauq * v * v**T and G(i) = I - taup * u * u**T */
/* > */
/* > where tauq and taup are real scalars, and v and u are real vectors;
*/
/* > v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
*/
/* > u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
*/
/* > tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* > The contents of A on exit are illustrated by the following examples: */
/* > */
/* > m = 6 and n = 5 (m > n): m = 5 and n = 6 (m < n): */
/* > */
/* > ( d e u1 u1 u1 ) ( d u1 u1 u1 u1 u1 ) */
/* > ( v1 d e u2 u2 ) ( e d u2 u2 u2 u2 ) */
/* > ( v1 v2 d e u3 ) ( v1 e d u3 u3 u3 ) */
/* > ( v1 v2 v3 d e ) ( v1 v2 e d u4 u4 ) */
/* > ( v1 v2 v3 v4 d ) ( v1 v2 v3 e d u5 ) */
/* > ( v1 v2 v3 v4 v5 ) */
/* > */
/* > where d and e denote diagonal and off-diagonal elements of B, vi */
/* > denotes an element of the vector defining H(i), and ui an element of */
/* > the vector defining G(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int lapack_sgebrd(integer *m, integer *n, real *a, integer *lda, real *d__, real *e, real *tauq, real *taup, real *work, integer * lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer i__, nbmin;
    extern /* Subroutine */
    void sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    integer minmn;
    extern /* Subroutine */
    int lapack_sgebd2(integer *, integer *, real *, integer *, real *, real *, real *, real *, real *, integer *);
    integer nb, nx;
    extern /* Subroutine */
    void slabrd_(integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, real *, integer *, real *, integer *);
    integer ws;
    extern /* Subroutine */
    void xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;
    /* Function Body */
    *info = 0;
    /* Computing MAX */
    i__1 = 1;
    i__2 = ilaenv_(&c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
    nb = fla_max(i__1,i__2);
    lwkopt = (*m + *n) * nb;
    work[1] = sroundup_lwork(&lwkopt);
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = fla_max(1,*m);
        if (*lwork < fla_max(i__1,*n) && ! lquery)
        {
            *info = -10;
        }
    }
    if (*info < 0)
    {
        i__1 = -(*info);
        xerbla_("SGEBRD", &i__1, (ftnlen)6);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    minmn = fla_min(*m,*n);
    if (minmn == 0)
    {
        work[1] = 1.f;
        return 0;
    }
    ws = fla_max(*m,*n);
    if (nb > 1 && nb < minmn)
    {
        /* Set the crossover point NX. */
        /* Computing MAX */
        i__1 = nb;
        i__2 = ilaenv_(&c__3, "SGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
        nx = fla_max(i__1,i__2);
        /* Determine when to switch from blocked to unblocked code. */
        if (nx < minmn)
        {
            ws = (*m + *n) * nb;
            if (*lwork < ws)
            {
                /* Not enough work space for the optimal NB, consider using */
                /* a smaller block size. */
                nbmin = ilaenv_(&c__2, "SGEBRD", " ", m, n, &c_n1, &c_n1);
                if (*lwork >= (*m + *n) * nbmin)
                {
                    nb = *lwork / (*m + *n);
                }
                else
                {
                    nb = 1;
                    nx = minmn;
                }
            }
        }
    }
    else
    {
        nx = minmn;
    }
    i__1 = minmn - nx;
    i__2 = nb;
    i__ = 1;

/* Current blocked algorithm has accuracy issue, so unblocked algorithm is enabled by default
Todo: This is a temporary workaround until the issue in the blocked algorithm is fixed.
*/
#if !FLA_ENABLE_AMD_OPT
    integer ldwrkx, ldwrky, j, i__4, i__3;
    static real c_b22 = 1.f;
    static real c_b21 = -1.f;
    ldwrkx = *m;
    ldwrky = *n;

    for (i__ = 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        /* Reduce rows and columns i:i+nb-1 to bidiagonal form and return */
        /* the matrices X and Y which are needed to update the unreduced */
        /* part of the matrix */
        i__3 = *m - i__ + 1;
        i__4 = *n - i__ + 1;
        slabrd_(&i__3, &i__4, &nb, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[ i__], &tauq[i__], &taup[i__], &work[1], &ldwrkx, &work[ldwrkx * nb + 1], &ldwrky);
        /* Update the trailing submatrix A(i+nb:m,i+nb:n), using an update */
        /* of the form A := A - V*Y**T - X*U**T */
        i__3 = *m - i__ - nb + 1;
        i__4 = *n - i__ - nb + 1;
        sgemm_("No transpose", "Transpose", &i__3, &i__4, &nb, &c_b21, &a[i__ + nb + i__ * a_dim1], lda, &work[ldwrkx * nb + nb + 1], & ldwrky, &c_b22, &a[i__ + nb + (i__ + nb) * a_dim1], lda);
        i__3 = *m - i__ - nb + 1;
        i__4 = *n - i__ - nb + 1;
        sgemm_("No transpose", "No transpose", &i__3, &i__4, &nb, &c_b21, & work[nb + 1], &ldwrkx, &a[i__ + (i__ + nb) * a_dim1], lda, & c_b22, &a[i__ + nb + (i__ + nb) * a_dim1], lda);
        /* Copy diagonal and off-diagonal elements of B back into A */
        if (*m >= *n)
        {
            i__3 = i__ + nb - 1;
            for (j = i__;
                    j <= i__3;
                    ++j)
            {
                a[j + j * a_dim1] = d__[j];
                a[j + (j + 1) * a_dim1] = e[j];
                /* L10: */
            }
        }
        else
        {
            i__3 = i__ + nb - 1;
            for (j = i__;
                    j <= i__3;
                    ++j)
            {
                a[j + j * a_dim1] = d__[j];
                a[j + 1 + j * a_dim1] = e[j];
                /* L20: */
            }
        }
        /* L30: */
    }
#else
    /* Use unblocked code to reduce the remainder of the matrix */
    i__2 = *m - i__ + 1;
    i__1 = *n - i__ + 1;
    lapack_sgebd2(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], & tauq[i__], &taup[i__], &work[1], info);
#endif
    work[1] = sroundup_lwork(&ws);
    return 0;
    /* End of SGEBRD */
}
/* lapack_sgebrd */
