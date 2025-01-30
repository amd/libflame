/* ./sgelq.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
/* > \brief \b SGELQ */
/* Definition: */
/* =========== */
/* SUBROUTINE SGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N, TSIZE, LWORK */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), T( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGELQ computes an LQ factorization of a real M-by-N matrix A: */
/* > */
/* > A = ( L 0 ) * Q */
/* > */
/* > where: */
/* > */
/* > Q is a N-by-N orthogonal matrix;
 */
/* > L is a lower-triangular M-by-M matrix;
 */
/* > 0 is a M-by-(N-M) zero matrix, if M < N. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the elements on and below the diagonal of the array */
/* > contain the M-by-min(M,N) lower trapezoidal matrix L */
/* > (L is lower triangular if M <= N);
 */
/* > the elements above the diagonal are used to store part of the */
/* > data structure to represent Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is REAL array, dimension (MAX(5,TSIZE)) */
/* > On exit, if INFO = 0, T(1) returns optimal (or either minimal */
/* > or optimal, if query is assumed) TSIZE. See TSIZE for details. */
/* > Remaining T contains part of the data structure used to represent Q. */
/* > If one wants to apply or construct Q, then one needs to keep T */
/* > (in addition to A) and pass it to further subroutines. */
/* > \endverbatim */
/* > */
/* > \param[in] TSIZE */
/* > \verbatim */
/* > TSIZE is INTEGER */
/* > If TSIZE >= 5, the dimension of the array T. */
/* > If TSIZE = -1 or -2, then a workspace query is assumed. The routine */
/* > only calculates the sizes of the T and WORK arrays, returns these */
/* > values as the first entries of the T and WORK arrays, and no error */
/* > message related to T or WORK is issued by XERBLA. */
/* > If TSIZE = -1, the routine calculates optimal size of T for the */
/* > optimum performance and returns this value in T(1). */
/* > If TSIZE = -2, the routine calculates minimal size of T and */
/* > returns this value in T(1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) contains optimal (or either minimal */
/* > or optimal, if query was assumed) LWORK. */
/* > See LWORK for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If LWORK = -1 or -2, then a workspace query is assumed. The routine */
/* > only calculates the sizes of the T and WORK arrays, returns these */
/* > values as the first entries of the T and WORK arrays, and no error */
/* > message related to T or WORK is issued by XERBLA. */
/* > If LWORK = -1, the routine calculates optimal size of WORK for the */
/* > optimal performance and returns this value in WORK(1). */
/* > If LWORK = -2, the routine calculates minimal size of WORK and */
/* > returns this value in WORK(1). */
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
/* > \par Further Details */
/* ==================== */
/* > */
/* > \verbatim */
/* > */
/* > The goal of the interface is to give maximum freedom to the developers for */
/* > creating any LQ factorization algorithm they wish. The triangular */
/* > (trapezoidal) L has to be stored in the lower part of A. The lower part of A */
/* > and the array T can be used to store any relevant information for applying or */
/* > constructing the Q factor. The WORK array can safely be discarded after exit. */
/* > */
/* > Caution: One should not expect the sizes of T and WORK to be the same from one */
/* > LAPACK implementation to the other, or even from one execution to the other. */
/* > A workspace query (for T and WORK) is needed at each execution. However, */
/* > for a given execution, the size of T and WORK are fixed and will not change */
/* > from one query to the next. */
/* > */
/* > \endverbatim */
/* > */
/* > \par Further Details particular to this LAPACK implementation: */
/* ============================================================== */
/* > */
/* > \verbatim */
/* > */
/* > These details are particular for this LAPACK implementation. Users should not */
/* > take them for granted. These details may change in the future, and are not likely */
/* > true for another LAPACK implementation. These details are relevant if one wants */
/* > to try to understand the code. They are not part of the interface. */
/* > */
/* > In this version, */
/* > */
/* > T(2): row block size (MB) */
/* > T(3): column block size (NB) */
/* > T(6:TSIZE): data structure needed for Q, computed by */
/* > SLASWLQ or SGELQT */
/* > */
/* > Depending on the matrix dimensions M and N, and row and column */
/* > block sizes MB and NB returned by ILAENV, SGELQ will use either */
/* > SLASWLQ (if the matrix is short-and-wide) or SGELQT to compute */
/* > the LQ factorization. */
/* > \endverbatim */
/* > */
/* > \ingroup gelq */
/* > */
/* ===================================================================== */
/* Subroutine */
void sgelq_(integer *m, integer *n, real *a, integer *lda, real *t, integer *tsize, real *work,
            integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgelq inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", tsize %" FLA_IS
                      "",
                      *m, *n, *lda, *tsize);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer mb, nb;
    logical mint, minw;
    integer lwmin, lwreq, lwopt, nblcks;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
        void
        sgelqt_(integer *, integer *, integer *, real *, integer *, real *, integer *, real *,
                integer *);
    logical lminws, lquery;
    integer mintsz;
    extern real sroundup_lwork(integer *);
    extern /* Subroutine */
        void
        slaswlq_(integer *, integer *, integer *, integer *, real *, integer *, real *, integer *,
                 real *, integer *, integer *);
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. -- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --t;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *tsize == -1 || *tsize == -2 || *lwork == -1 || *lwork == -2;
    mint = FALSE_;
    minw = FALSE_;
    if(*tsize == -2 || *lwork == -2)
    {
        if(*tsize != -1)
        {
            mint = TRUE_;
        }
        if(*lwork != -1)
        {
            minw = TRUE_;
        }
    }
    /* Determine the block size */
    if(fla_min(*m, *n) > 0)
    {
        mb = ilaenv_(&c__1, "SGELQ ", " ", m, n, &c__1, &c_n1);
        nb = ilaenv_(&c__1, "SGELQ ", " ", m, n, &c__2, &c_n1);
    }
    else
    {
        mb = 1;
        nb = *n;
    }
    if(mb > fla_min(*m, *n) || mb < 1)
    {
        mb = 1;
    }
    if(nb > *n || nb <= *m)
    {
        nb = *n;
    }
    mintsz = *m + 5;
    if(nb > *m && *n > *m)
    {
        if((*n - *m) % (nb - *m) == 0)
        {
            nblcks = (*n - *m) / (nb - *m);
        }
        else
        {
            nblcks = (*n - *m) / (nb - *m) + 1;
        }
    }
    else
    {
        nblcks = 1;
    }
    /* Determine if the workspace size satisfies minimal size */
    if(*n <= *m || nb <= *m || nb >= *n)
    {
        lwmin = fla_max(1, *n);
        /* Computing MAX */
        i__1 = 1;
        i__2 = mb * *n; // , expr subst
        lwopt = fla_max(i__1, i__2);
    }
    else
    {
        lwmin = fla_max(1, *m);
        /* Computing MAX */
        i__1 = 1;
        i__2 = mb * *m; // , expr subst
        lwopt = fla_max(i__1, i__2);
    }
    lminws = FALSE_;
    /* Computing MAX */
    i__1 = 1;
    i__2 = mb * *m * nblcks + 5; // , expr subst
    if((*tsize < fla_max(i__1, i__2) || *lwork < lwopt) && *lwork >= lwmin && *tsize >= mintsz
       && !lquery)
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = mb * *m * nblcks + 5; // , expr subst
        if(*tsize < fla_max(i__1, i__2))
        {
            lminws = TRUE_;
            mb = 1;
            nb = *n;
        }
        if(*lwork < lwopt)
        {
            lminws = TRUE_;
            mb = 1;
        }
    }
    if(*n <= *m || nb <= *m || nb >= *n)
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = mb * *n; // , expr subst
        lwreq = fla_max(i__1, i__2);
    }
    else
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = mb * *m; // , expr subst
        lwreq = fla_max(i__1, i__2);
    }
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = mb * *m * nblcks + 5; // , expr subst
        if(*tsize < fla_max(i__1, i__2) && !lquery && !lminws)
        {
            *info = -6;
        }
        else if(*lwork < lwreq && !lquery && !lminws)
        {
            *info = -8;
        }
    }
    if(*info == 0)
    {
        if(mint)
        {
            t[1] = (real)mintsz;
        }
        else
        {
            t[1] = (real)(mb * *m * nblcks + 5);
        }
        t[2] = (real)mb;
        t[3] = (real)nb;
        if(minw)
        {
            work[1] = sroundup_lwork(&lwmin);
        }
        else
        {
            work[1] = sroundup_lwork(&lwreq);
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGELQ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(fla_min(*m, *n) == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* The LQ Decomposition */
    if(*n <= *m || nb <= *m || nb >= *n)
    {
        sgelqt_(m, n, &mb, &a[a_offset], lda, &t[6], &mb, &work[1], info);
    }
    else
    {
        slaswlq_(m, n, &mb, &nb, &a[a_offset], lda, &t[6], &mb, &work[1], lwork, info);
    }
    work[1] = sroundup_lwork(&lwreq);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGELQ */
}
/* sgelq_ */
