/* ../netlib/dgetri.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link
 with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in
 that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
#include "fla_lapack_x86_common.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b20 = -1.;
static doublereal c_b22 = 1.;
/* > \brief \b DGETRI */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGETRI + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetri.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetri.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetri.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGETRI computes the inverse of a matrix using the LU factorization */
/* > computed by DGETRF. */
/* > */
/* > This method inverts U and then computes inv(A) by solving the system */
/* > inv(A)*L = inv(U) for inv(A). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the factors L and U from the factorization */
/* > A = P*L*U as computed by DGETRF. */
/* > On exit, if INFO = 0, the inverse of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from DGETRF;
for 1<=i<=N, row i of the */
/* > matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO=0, then WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
/* > For optimal performance LWORK >= N*NB, where NB is */
/* > the optimal blocksize returned by ILAENV. */
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
/* > > 0: if INFO = i, U(i,i) is exactly zero;
the matrix is */
/* > singular and its inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup doubleGEcomputational */
/* ===================================================================== */
/* Subroutine */
void dgetri_(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work,
             integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgetri inputs: n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS "", *n, *lda,
                      *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j, jb, nb = 0, jj, jp, nn, iws;
    extern /* Subroutine */
        void
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *),
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *);
    integer nbmin;
    extern /* Subroutine */
        void
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *),
        dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ldwork;
    extern /* Subroutine */
        void
        dtrtri_(char *, char *, integer *, doublereal *, integer *, integer *);
    integer lwkopt;
    logical lquery;
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
    --ipiv;
    --work;
    /* Function Body */
    *info = 0;
    if(*n > 0 && *n <= 16)
    {
        lwkopt = *n;
    }
    else
    {
        nb = ilaenv_(&c__1, "DGETRI", " ", n, &c_n1, &c_n1, &c_n1);
        lwkopt = *n * nb;
    }
    work[1] = (doublereal)lwkopt;
    lquery = *lwork == -1;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -3;
    }
    else if(*lwork < fla_max(1, *n) && !lquery)
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGETRI", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Form inv(U). If INFO > 0 from DTRTRI, then U is singular, */
    /* and the inverse is not computed. */
#if FLA_ENABLE_AMD_OPT
    if(*n > 0 && *n <= 16)
    {
        lapack_getri_small_d(n, a, lda, ipiv, work, info);
    }
    else
#endif
    {
        dtrtri_("Upper", "Non-unit", n, &a[a_offset], lda, info);

        if(*info > 0)
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        nbmin = 2;
        ldwork = *n;
        if(nb > 1 && nb < *n)
        {
            /* Computing MAX */
            i__1 = ldwork * nb;
            iws = fla_max(i__1, 1);
            if(*lwork < iws)
            {
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = ilaenv_(&c__2, "DGETRI", " ", n, &c_n1, &c_n1, &c_n1); // , expr subst
                nbmin = fla_max(i__1, i__2);
            }
        }
        else
        {
            iws = *n;
        }
        /* Solve the equation inv(A)*L = inv(U) for inv(A). */
        if(nb < nbmin || nb >= *n)
        {
            /* Use unblocked code. */
            for(j = *n; j >= 1; --j)
            {
                /* Copy current column of L to WORK and replace with zeros. */
                i__1 = *n;
                for(i__ = j + 1; i__ <= i__1; ++i__)
                {
                    work[i__] = a[i__ + j * a_dim1];
                    a[i__ + j * a_dim1] = 0.;
                    /* L10: */
                }
                /* Compute current column of inv(A). */
                if(j < *n)
                {
                    i__1 = *n - j;
                    dgemv_("No transpose", n, &i__1, &c_b20, &a[(j + 1) * a_dim1 + 1], lda,
                           &work[j + 1], &c__1, &c_b22, &a[j * a_dim1 + 1], &c__1);
                }
                /* L20: */
            }
        }
        else
        {
            /* Use blocked code. */
            nn = (*n - 1) / nb * nb + 1;
            i__1 = -nb;
            for(j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1)
            {
                /* Computing MIN */
                i__2 = nb;
                i__3 = *n - j + 1; // , expr subst
                jb = fla_min(i__2, i__3);
                /* Copy current block column of L to WORK and replace with */
                /* zeros. */
                i__2 = j + jb - 1;
                for(jj = j; jj <= i__2; ++jj)
                {
                    i__3 = *n;
                    for(i__ = jj + 1; i__ <= i__3; ++i__)
                    {
                        work[i__ + (jj - j) * ldwork] = a[i__ + jj * a_dim1];
                        a[i__ + jj * a_dim1] = 0.;
                        /* L30: */
                    }
                    /* L40: */
                }
                /* Compute current block column of inv(A). */
                if(j + jb <= *n)
                {
                    i__2 = *n - j - jb + 1;
                    dgemm_("No transpose", "No transpose", n, &jb, &i__2, &c_b20,
                           &a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &ldwork, &c_b22,
                           &a[j * a_dim1 + 1], lda);
                }
                dtrsm_("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b22, &work[j], &ldwork,
                       &a[j * a_dim1 + 1], lda);
                /* L50: */
            }
        }
        /* Apply column interchanges. */
        for(j = *n - 1; j >= 1; --j)
        {
            jp = ipiv[j];
            if(jp != j)
            {
                dswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
            }
            /* L60: */
        }
        work[1] = (doublereal)iws;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DGETRI */
}
/* dgetri_ */
