/* ../netlib/zhecon_rook.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief <b> ZHECON_ROOK estimates the reciprocal of the condition number fort HE matrices using
 * factorizat ion obtained with one of the bounded diagonal pivoting methods (max 2 interchanges)
 * </b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHECON_ROOK + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhecon_
 * rook.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhecon_
 * rook.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhecon_
 * rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHECON_ROOK( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* DOUBLE PRECISION ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHECON_ROOK estimates the reciprocal of the condition number of a complex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by CHETRF_ROOK. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**H;
 */
/* > = 'L': Lower triangular, form is A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by CHETRF_ROOK. */
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
/* > Details of the interchanges and the block structure of D */
/* > as determined by CHETRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is DOUBLE PRECISION */
/* > The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is DOUBLE PRECISION */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* > estimate of the 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (2*N) */
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
/* > \date June 2017 */
/* > \ingroup complex16HEcomputational */
/* > \par Contributors: */
/* ================== */
/* > \verbatim */
/* > */
/* > June 2017, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* > School of Mathematics, */
/* > University of Manchester */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
void zhecon_rook_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv,
                  doublereal *anorm, doublereal *rcond, doublecomplex *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhecon_rook inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n,
                      *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer i__;
    extern /* Subroutine */
        void
        zhetrs_rook_(char *, integer *, integer *, doublecomplex *, integer *, integer *,
                     doublecomplex *, integer *, integer *);
    integer kase;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    logical upper;
    extern /* Subroutine */
        void
        zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal ainvnm;
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
    /* .. Local Arrays .. */
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
    else if(*anorm < 0.)
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHECON_ROOK", &i__1, (ftnlen)11);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *rcond = 0.;
    if(*n == 0)
    {
        *rcond = 1.;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*anorm <= 0.)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if(upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for(i__ = *n; i__ >= 1; --i__)
        {
            i__1 = i__ + i__ * a_dim1;
            if(ipiv[i__] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* L10: */
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = i__ + i__ * a_dim1;
            if(ipiv[i__] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* L20: */
        }
    }
    /* Estimate the 1-norm of the inverse. */
    kase = 0;
L30:
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if(kase != 0)
    {
        /* Multiply by inv(L*D*L**H) or inv(U*D*U**H). */
        zhetrs_rook_(uplo, n, &c__1, &a[a_offset], lda, &ipiv[1], &work[1], n, info);
        goto L30;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if(ainvnm != 0.)
    {
        *rcond = 1. / ainvnm / *anorm;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZHECON_ROOK */
}
/* zhecon_rook__ */
