/* ./zgebal.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZGEBAL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGEBAL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebal.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebal.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebal.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOB */
/* INTEGER IHI, ILO, INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION SCALE( * ) */
/* COMPLEX*16 A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEBAL balances a general complex matrix A. This involves, first, */
/* > permuting A by a similarity transformation to isolate eigenvalues */
/* > in the first 1 to ILO-1 and last IHI+1 to N elements on the */
/* > diagonal;
and second, applying a diagonal similarity transformation */
/* > to rows and columns ILO to IHI to make the rows and columns as */
/* > close in norm as possible. Both steps are optional. */
/* > */
/* > Balancing may reduce the 1-norm of the matrix, and improve the */
/* > accuracy of the computed eigenvalues and/or eigenvectors. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is CHARACTER*1 */
/* > Specifies the operations to be performed on A: */
/* > = 'N': none: simply set ILO = 1, IHI = N, SCALE(I) = 1.0 */
/* > for i = 1,...,N;
 */
/* > = 'P': permute only;
 */
/* > = 'S': scale only;
 */
/* > = 'B': both permute and scale. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the input matrix A. */
/* > On exit, A is overwritten by the balanced matrix. */
/* > If JOB = 'N', A is not referenced. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > ILO and IHI are set to integers such that on exit */
/* > A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N. */
/* > If JOB = 'N' or 'S', ILO = 1 and IHI = N. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION array, dimension (N) */
/* > Details of the permutations and scaling factors applied to */
/* > A. If P(j) is the index of the row and column interchanged */
/* > with row and column j and D(j) is the scaling factor */
/* > applied to row and column j, then */
/* > SCALE(j) = P(j) for j = 1,...,ILO-1 */
/* > = D(j) for j = ILO,...,IHI */
/* > = P(j) for j = IHI+1,...,N. */
/* > The order in which the interchanges are made is N to IHI+1, */
/* > then 1 to ILO-1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gebal */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The permutations consist of row and column interchanges which put */
/* > the matrix in the form */
/* > */
/* > ( T1 X Y ) */
/* > P A P = ( 0 B Z ) */
/* > ( 0 0 T2 ) */
/* > */
/* > where T1 and T2 are upper triangular matrices whose eigenvalues lie */
/* > along the diagonal. The column indices ILO and IHI mark the starting */
/* > and ending columns of the submatrix B. Balancing consists of applying */
/* > a diagonal similarity transformation inv(D) * B * D to make the */
/* > 1-norms of each row of B and its corresponding column nearly equal. */
/* > The output matrix is */
/* > */
/* > ( T1 X*D Y ) */
/* > ( 0 inv(D)*B*D inv(D)*Z ). */
/* > ( 0 0 T2 ) */
/* > */
/* > Information about the permutations P and the diagonal matrix D is */
/* > returned in the vector SCALE. */
/* > */
/* > This subroutine is based on the EISPACK routine CBAL. */
/* > */
/* > Modified by Tzu-Yi Chen, Computer Science Division, University of */
/* > California at Berkeley, USA */
/* > */
/* > Refactored by Evert Provoost, Department of Computer Science, */
/* > KU Leuven, Belgium */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void zgebal_(char *job, integer *n, doublecomplex *a, integer *lda, integer *ilo, integer *ihi,
             doublereal *scale, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgebal inputs: job %c, n %" FLA_IS ", lda %" FLA_IS "", *job, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *);
    /* Local variables */
    doublereal c__, f, g;
    integer i__, j, k, l;
    doublereal r__, s, ca, ra;
    integer ica, ira;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(char *);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    logical noconv;
    logical canswap;
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
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --scale;
    /* Function Body */
    *info = 0;
    if(!lsame_(job, "N", 1, 1) && !lsame_(job, "P", 1, 1) && !lsame_(job, "S", 1, 1)
       && !lsame_(job, "B", 1, 1))
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
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGEBAL", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick returns. */
    if(*n == 0)
    {
        *ilo = 1;
        *ihi = 0;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(lsame_(job, "N", 1, 1))
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            scale[i__] = 1.;
        }
        *ilo = 1;
        *ihi = *n;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Permutation to isolate eigenvalues if possible. */
    k = 1;
    l = *n;
    if(!lsame_(job, "S", 1, 1))
    {
        /* Row and column exchange. */
        noconv = TRUE_;
        while(noconv)
        {
            /* Search for rows isolating an eigenvalue and push them down. */
            noconv = FALSE_;
            for(i__ = l; i__ >= 1; --i__)
            {
                canswap = TRUE_;
                i__1 = l;
                for(j = 1; j <= i__1; ++j)
                {
                    i__2 = i__ + j * a_dim1;
                    if(i__ != j && (a[i__2].r != 0. || d_imag(&a[i__ + j * a_dim1]) != 0.))
                    {
                        canswap = FALSE_;
                        break;
                    }
                }
                if(canswap)
                {
                    scale[l] = (doublereal)i__;
                    if(i__ != l)
                    {
                        zswap_(&l, &a[i__ * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
                        i__1 = *n - k + 1;
                        zswap_(&i__1, &a[i__ + k * a_dim1], lda, &a[l + k * a_dim1], lda);
                    }
                    noconv = TRUE_;
                    if(l == 1)
                    {
                        *ilo = 1;
                        *ihi = 1;
                        AOCL_DTL_TRACE_LOG_EXIT
                        return;
                    }
                    --l;
                }
            }
        }
        noconv = TRUE_;
        while(noconv)
        {
            /* Search for columns isolating an eigenvalue and push them left. */
            noconv = FALSE_;
            i__1 = l;
            for(j = k; j <= i__1; ++j)
            {
                canswap = TRUE_;
                i__2 = l;
                for(i__ = k; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * a_dim1;
                    if(i__ != j && (a[i__3].r != 0. || d_imag(&a[i__ + j * a_dim1]) != 0.))
                    {
                        canswap = FALSE_;
                        break;
                    }
                }
                if(canswap)
                {
                    scale[k] = (doublereal)j;
                    if(j != k)
                    {
                        zswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
                        i__2 = *n - k + 1;
                        zswap_(&i__2, &a[j + k * a_dim1], lda, &a[k + k * a_dim1], lda);
                    }
                    noconv = TRUE_;
                    ++k;
                }
            }
        }
    }
    /* Initialize SCALE for non-permuted submatrix. */
    i__1 = l;
    for(i__ = k; i__ <= i__1; ++i__)
    {
        scale[i__] = 1.;
    }
    /* If we only had to permute, we are done. */
    if(lsame_(job, "P", 1, 1))
    {
        *ilo = k;
        *ihi = l;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Balance the submatrix in rows K to L. */
    /* Iterative loop for norm reduction. */
    sfmin1 = dlamch_("S") / dlamch_("P");
    sfmax1 = 1. / sfmin1;
    sfmin2 = sfmin1 * 2.;
    sfmax2 = 1. / sfmin2;
    noconv = TRUE_;
    while(noconv)
    {
        noconv = FALSE_;
        i__1 = l;
        for(i__ = k; i__ <= i__1; ++i__)
        {
            i__2 = l - k + 1;
            c__ = dznrm2_(&i__2, &a[k + i__ * a_dim1], &c__1);
            i__2 = l - k + 1;
            r__ = dznrm2_(&i__2, &a[i__ + k * a_dim1], lda);
            ica = izamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
            ca = z_abs(&a[ica + i__ * a_dim1]);
            i__2 = *n - k + 1;
            ira = izamax_(&i__2, &a[i__ + k * a_dim1], lda);
            ra = z_abs(&a[i__ + (ira + k - 1) * a_dim1]);
            /* Guard against zero C or R due to underflow. */
            if(c__ == 0. || r__ == 0.)
            {
                continue;
            }
            /* Exit if NaN to avoid infinite loop */
            d__1 = c__ + ca + r__ + ra;
            if(disnan_(&d__1))
            {
                *info = -3;
                i__2 = -(*info);
                xerbla_("ZGEBAL", &i__2, (ftnlen)6);
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            g = r__ / 2.;
            f = 1.;
            s = c__ + r__;
            for(;;)
            {
                /* while(complicated condition) */
                /* Computing MAX */
                d__1 = fla_max(f, c__);
                /* Computing MIN */
                d__2 = fla_min(r__, g);
                if(!(c__ < g && fla_max(d__1, ca) < sfmax2 && fla_min(d__2, ra) > sfmin2))
                    break;
                f *= 2.;
                c__ *= 2.;
                ca *= 2.;
                r__ /= 2.;
                g /= 2.;
                ra /= 2.;
            }
            g = c__ / 2.;
            for(;;)
            {
                /* while(complicated condition) */
                /* Computing MIN */
                d__1 = fla_min(f, c__);
                d__1 = fla_min(d__1, g); // , expr subst
                if(!(g >= r__ && fla_max(r__, ra) < sfmax2 && fla_min(d__1, ca) > sfmin2))
                    break;
                f /= 2.;
                c__ /= 2.;
                g /= 2.;
                ca /= 2.;
                r__ *= 2.;
                ra *= 2.;
            }
            /* Now balance. */
            if(c__ + r__ >= s * .95)
            {
                continue;
            }
            if(f < 1. && scale[i__] < 1.)
            {
                if(f * scale[i__] <= sfmin1)
                {
                    continue;
                }
            }
            if(f > 1. && scale[i__] > 1.)
            {
                if(scale[i__] >= sfmax1 / f)
                {
                    continue;
                }
            }
            g = 1. / f;
            scale[i__] *= f;
            noconv = TRUE_;
            i__2 = *n - k + 1;
            zdscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
            zdscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);
        }
    }
    *ilo = k;
    *ihi = l;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGEBAL */
}
/* zgebal_ */
