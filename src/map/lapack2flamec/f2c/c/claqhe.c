/* ../netlib/claqhe.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLAQHE scales a Hermitian matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQHE + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqhe.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqhe.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqhe.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, UPLO */
/* INTEGER LDA, N */
/* REAL AMAX, SCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL S( * ) */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQHE equilibrates a Hermitian matrix A using the scaling factors */
/* > in the vector S. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > Hermitian matrix A is stored. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > n by n upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading n by n lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if EQUED = 'Y', the equilibrated matrix: */
/* > diag(S) * A * diag(S). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is REAL array, dimension (N) */
/* > The scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] SCOND */
/* > \verbatim */
/* > SCOND is REAL */
/* > Ratio of the smallest S(i) to the largest S(i). */
/* > \endverbatim */
/* > */
/* > \param[in] AMAX */
/* > \verbatim */
/* > AMAX is REAL */
/* > Absolute value of largest matrix entry. */
/* > \endverbatim */
/* > */
/* > \param[out] EQUED */
/* > \verbatim */
/* > EQUED is CHARACTER*1 */
/* > Specifies whether or not equilibration was done. */
/* > = 'N': No equilibration. */
/* > = 'Y': Equilibration was done, i.e., A has been replaced by */
/* > diag(S) * A * diag(S). */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > THRESH is a threshold value used to decide if scaling should be done */
/* > based on the ratio of the scaling factors. If SCOND < THRESH, */
/* > scaling is done. */
/* > */
/* > LARGE and SMALL are threshold values used to decide if scaling should */
/* > be done based on the absolute size of the largest matrix element. */
/* > If AMAX > LARGE or AMAX < SMALL, scaling is done. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexHEauxiliary */
/* ===================================================================== */
/* Subroutine */
void claqhe_(char *uplo, integer *n, complex *a, integer *lda, real *s, real *scond, real *amax,
             char *equed)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "claqhe inputs: uplo %c, n %lld, lda %lld", *uplo, *n, *lda);
#else
    snprintf(buffer, 256, "claqhe inputs: uplo %c, n %d, lda %d", *uplo, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1;
    complex q__1;
    /* Local variables */
    integer i__, j;
    real cj, large;
    extern logical lsame_(char *, char *, integer, integer);
    real small_val;
    extern real slamch_(char *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    /* Function Body */
    if(*n <= 0)
    {
        *(unsigned char *)equed = 'N';
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Initialize LARGE and SMALL. */
    small_val = slamch_("Safe minimum") / slamch_("Precision");
    large = 1.f / small_val;
    if(*scond >= .1f && *amax >= small_val && *amax <= large)
    {
        /* No equilibration */
        *(unsigned char *)equed = 'N';
    }
    else
    {
        /* Replace A by diag(S) * A * diag(S). */
        if(lsame_(uplo, "U", 1, 1))
        {
            /* Upper triangle of A is stored. */
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                cj = s[j];
                i__2 = j - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * a_dim1;
                    r__1 = cj * s[i__];
                    i__4 = i__ + j * a_dim1;
                    q__1.r = r__1 * a[i__4].r;
                    q__1.i = r__1 * a[i__4].i; // , expr subst
                    a[i__3].r = q__1.r;
                    a[i__3].i = q__1.i; // , expr subst
                    /* L10: */
                }
                i__2 = j + j * a_dim1;
                i__3 = j + j * a_dim1;
                r__1 = cj * cj * a[i__3].r;
                a[i__2].r = r__1;
                a[i__2].i = 0.f; // , expr subst
                /* L20: */
            }
        }
        else
        {
            /* Lower triangle of A is stored. */
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                cj = s[j];
                i__2 = j + j * a_dim1;
                i__3 = j + j * a_dim1;
                r__1 = cj * cj * a[i__3].r;
                a[i__2].r = r__1;
                a[i__2].i = 0.f; // , expr subst
                i__2 = *n;
                for(i__ = j + 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * a_dim1;
                    r__1 = cj * s[i__];
                    i__4 = i__ + j * a_dim1;
                    q__1.r = r__1 * a[i__4].r;
                    q__1.i = r__1 * a[i__4].i; // , expr subst
                    a[i__3].r = q__1.r;
                    a[i__3].i = q__1.i; // , expr subst
                    /* L30: */
                }
                /* L40: */
            }
        }
        *(unsigned char *)equed = 'Y';
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLAQHE */
}
/* claqhe_ */
