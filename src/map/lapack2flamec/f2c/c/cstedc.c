/* ./cstedc.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static real c_b17 = 0.f;
static real c_b18 = 1.f;
static integer c__1 = 1;
/* > \brief \b CSTEDC */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSTEDC + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cstedc.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cstedc.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cstedc.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, */
/* LRWORK, IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER COMPZ */
/* INTEGER INFO, LDZ, LIWORK, LRWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL D( * ), E( * ), RWORK( * ) */
/* COMPLEX WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSTEDC computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
/* > The eigenvectors of a full or band complex Hermitian matrix can also */
/* > be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this */
/* > matrix to tridiagonal form. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] COMPZ */
/* > \verbatim */
/* > COMPZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only. */
/* > = 'I': Compute eigenvectors of tridiagonal matrix also. */
/* > = 'V': Compute eigenvectors of original Hermitian matrix */
/* > also. On entry, Z contains the unitary matrix used */
/* > to reduce the original matrix to tridiagonal form. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the diagonal elements of the tridiagonal matrix. */
/* > On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > On entry, the subdiagonal elements of the tridiagonal matrix. */
/* > On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ,N) */
/* > On entry, if COMPZ = 'V', then Z contains the unitary */
/* > matrix used in the reduction to tridiagonal form. */
/* > On exit, if INFO = 0, then if COMPZ = 'V', Z contains the */
/* > orthonormal eigenvectors of the original Hermitian matrix, */
/* > and if COMPZ = 'I', Z contains the orthonormal eigenvectors */
/* > of the symmetric tridiagonal matrix. */
/* > If COMPZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1. */
/* > If eigenvectors are desired, then LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1. */
/* > If COMPZ = 'V' and N > 1, LWORK must be at least N*N. */
/* > Note that for COMPZ = 'V', then if N is less than or */
/* > equal to the minimum divide size, usually 25, then LWORK need */
/* > only be 1. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal sizes of the WORK, RWORK and */
/* > IWORK arrays, returns these values as the first entries of */
/* > the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* > On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of the array RWORK. */
/* > If COMPZ = 'N' or N <= 1, LRWORK must be at least 1. */
/* > If COMPZ = 'V' and N > 1, LRWORK must be at least */
/* > 1 + 3*N + 2*N*lg N + 4*N**2 , */
/* > where lg( N ) = smallest integer k such */
/* > that 2**k >= N. */
/* > If COMPZ = 'I' and N > 1, LRWORK must be at least */
/* > 1 + 4*N + 2*N**2 . */
/* > Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* > equal to the minimum divide size, usually 25, then LRWORK */
/* > need only be fla_max(1,2*(N-1)). */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of the array IWORK. */
/* > If COMPZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* > If COMPZ = 'V' or N > 1, LIWORK must be at least */
/* > 6 + 6*N + 5*N*lg N. */
/* > If COMPZ = 'I' or N > 1, LIWORK must be at least */
/* > 3 + 5*N . */
/* > Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* > equal to the minimum divide size, usually 25, then LIWORK */
/* > need only be 1. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: The algorithm failed to compute an eigenvalue while */
/* > working on the submatrix lying in rows and columns */
/* > INFO/(N+1) through mod(INFO,N+1). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup stedc */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* ===================================================================== */
/* Subroutine */
void cstedc_(char *compz, integer *n, real *d__, real *e, complex *z__, integer *ldz, complex *work,
             integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *liwork,
             integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cstedc inputs: compz %c, n %lld, ldz %lld, lwork %lld, lrwork %lld, liwork %lld",
             *compz, *n, *ldz, *lwork, *lrwork, *liwork);
#else
    snprintf(buffer, 256, "cstedc inputs: compz %c, n %d, ldz %d, lwork %d, lrwork %d, liwork %d",
             *compz, *n, *ldz, *lwork, *lrwork, *liwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, m;
    real p;
    integer ii, ll, lgn;
    real eps, tiny;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        cswap_(integer *, complex *, integer *, complex *, integer *);
    integer lwmin;
    extern /* Subroutine */
        void
        claed0_(integer *, integer *, real *, real *, complex *, integer *, complex *, integer *,
                real *, integer *, integer *);
    integer start;
    extern /* Subroutine */
        void
        clacrm_(integer *, integer *, complex *, integer *, real *, integer *, complex *, integer *,
                real *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer finish;
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *),
        sstedc_(char *, integer *, real *, real *, real *, integer *, real *, integer *, integer *,
                integer *, integer *),
        slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    integer liwmin, icompz;
    extern /* Subroutine */
        void
        csteqr_(char *, integer *, real *, real *, complex *, integer *, real *, integer *);
    real orgnrm;
    extern real slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */
        void
        ssterf_(integer *, real *, real *, integer *);
    integer lrwmin;
    logical lquery;
    integer smlsiz;
    extern /* Subroutine */
        void
        ssteqr_(char *, integer *, real *, real *, real *, integer *, real *, integer *);
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    *info = 0;
    lwmin = 0;
    lrwmin = 0;
    liwmin = 0;
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;
    if(lsame_(compz, "N", 1, 1))
    {
        icompz = 0;
    }
    else if(lsame_(compz, "V", 1, 1))
    {
        icompz = 1;
    }
    else if(lsame_(compz, "I", 1, 1))
    {
        icompz = 2;
    }
    else
    {
        icompz = -1;
    }
    if(icompz < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*ldz < 1 || icompz > 0 && *ldz < fla_max(1, *n))
    {
        *info = -6;
    }
    if(*info == 0)
    {
        /* Compute the workspace requirements */
        smlsiz = ilaenv_(&c__9, "CSTEDC", " ", &c__0, &c__0, &c__0, &c__0);
        if(*n <= 1 || icompz == 0)
        {
            lwmin = 1;
            liwmin = 1;
            lrwmin = 1;
        }
        else if(*n <= smlsiz)
        {
            lwmin = 1;
            liwmin = 1;
            lrwmin = *n - 1 << 1;
        }
        else if(icompz == 1)
        {
            lgn = (integer)(log((real)(*n)) / log(2.f));
            if(pow_ii(&c__2, &lgn) < *n)
            {
                ++lgn;
            }
            if(pow_ii(&c__2, &lgn) < *n)
            {
                ++lgn;
            }
            lwmin = *n * *n;
            /* Computing 2nd power */
            i__1 = *n;
            lrwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
            liwmin = *n * 6 + 6 + *n * 5 * lgn;
        }
        else if(icompz == 2)
        {
            lwmin = 1;
            /* Computing 2nd power */
            i__1 = *n;
            lrwmin = (*n << 2) + 1 + (i__1 * i__1 << 1);
            liwmin = *n * 5 + 3;
        }
        r__1 = sroundup_lwork(&lwmin);
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
        rwork[1] = (real)lrwmin;
        iwork[1] = liwmin;
        if(*lwork < lwmin && !lquery)
        {
            *info = -8;
        }
        else if(*lrwork < lrwmin && !lquery)
        {
            *info = -10;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -12;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CSTEDC", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*n == 1)
    {
        if(icompz != 0)
        {
            i__1 = z_dim1 + 1;
            z__[i__1].r = 1.f;
            z__[i__1].i = 0.f; // , expr subst
        }
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* If the following conditional clause is removed, then the routine */
    /* will use the Divide and Conquer routine to compute only the */
    /* eigenvalues, which requires (3N + 3N**2) real workspace and */
    /* (2 + 5N + 2N lg(N)) integer workspace. */
    /* Since on many architectures SSTERF is much faster than any other */
    /* algorithm for finding eigenvalues only, it is used here */
    /* as the default. If the conditional clause is removed, then */
    /* information on the size of workspace needs to be changed. */
    /* If COMPZ = 'N', use SSTERF to compute the eigenvalues. */
    if(icompz == 0)
    {
        ssterf_(n, &d__[1], &e[1], info);
        goto L70;
    }
    /* If N is smaller than the minimum divide size (SMLSIZ+1), then */
    /* solve the problem with another solver. */
    if(*n <= smlsiz)
    {
        csteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &rwork[1], info);
    }
    else
    {
        /* If COMPZ = 'I', we simply call SSTEDC instead. */
        if(icompz == 2)
        {
            slaset_("Full", n, n, &c_b17, &c_b18, &rwork[1], n);
            ll = *n * *n + 1;
            i__1 = *lrwork - ll + 1;
            sstedc_("I", n, &d__[1], &e[1], &rwork[1], n, &rwork[ll], &i__1, &iwork[1], liwork,
                    info);
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * z_dim1;
                    i__4 = (j - 1) * *n + i__;
                    z__[i__3].r = rwork[i__4];
                    z__[i__3].i = 0.f; // , expr subst
                    /* L10: */
                }
                /* L20: */
            }
            goto L70;
        }
        /* From now on, only option left to be handled is COMPZ = 'V', */
        /* i.e. ICOMPZ = 1. */
        /* Scale. */
        orgnrm = slanst_("M", n, &d__[1], &e[1]);
        if(orgnrm == 0.f)
        {
            goto L70;
        }
        eps = slamch_("Epsilon");
        start = 1;
    /* while ( START <= N ) */
    L30:
        if(start <= *n)
        {
            /* Let FINISH be the position of the next subdiagonal entry */
            /* such that E( FINISH ) <= TINY or FINISH = N if no such */
            /* subdiagonal exists. The matrix identified by the elements */
            /* between START and FINISH constitutes an independent */
            /* sub-problem. */
            finish = start;
        L40:
            if(finish < *n)
            {
                tiny = eps * sqrt((r__1 = d__[finish], f2c_abs(r__1)))
                       * sqrt((r__2 = d__[finish + 1], f2c_abs(r__2)));
                if((r__1 = e[finish], f2c_abs(r__1)) > tiny)
                {
                    ++finish;
                    goto L40;
                }
            }
            /* (Sub) Problem determined. Compute its size and solve it. */
            m = finish - start + 1;
            if(m > smlsiz)
            {
                /* Scale. */
                orgnrm = slanst_("M", &m, &d__[start], &e[start]);
                slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[start], &m, info);
                i__1 = m - 1;
                i__2 = m - 1;
                slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[start], &i__2, info);
                claed0_(n, &m, &d__[start], &e[start], &z__[start * z_dim1 + 1], ldz, &work[1], n,
                        &rwork[1], &iwork[1], info);
                if(*info > 0)
                {
                    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m + 1) + start - 1;
                    goto L70;
                }
                /* Scale back. */
                slascl_("G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[start], &m, info);
            }
            else
            {
                ssteqr_("I", &m, &d__[start], &e[start], &rwork[1], &m, &rwork[m * m + 1], info);
                clacrm_(n, &m, &z__[start * z_dim1 + 1], ldz, &rwork[1], &m, &work[1], n,
                        &rwork[m * m + 1]);
                clacpy_("A", n, &m, &work[1], n, &z__[start * z_dim1 + 1], ldz);
                if(*info > 0)
                {
                    *info = start * (*n + 1) + finish;
                    goto L70;
                }
            }
            start = finish + 1;
            goto L30;
        }
        /* endwhile */
        /* Use Selection Sort to minimize swaps of eigenvectors */
        i__1 = *n;
        for(ii = 2; ii <= i__1; ++ii)
        {
            i__ = ii - 1;
            k = i__;
            p = d__[i__];
            i__2 = *n;
            for(j = ii; j <= i__2; ++j)
            {
                if(d__[j] < p)
                {
                    k = j;
                    p = d__[j];
                }
                /* L50: */
            }
            if(k != i__)
            {
                d__[k] = d__[i__];
                d__[i__] = p;
                cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1], &c__1);
            }
            /* L60: */
        }
    }
L70:
    r__1 = sroundup_lwork(&lwmin);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    rwork[1] = (real)lrwmin;
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CSTEDC */
}
/* cstedc_ */
