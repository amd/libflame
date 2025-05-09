/* ../netlib/zlaqsp.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLAQSP scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by sppequ. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAQSP + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqsp.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqsp.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqsp.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, UPLO */
/* INTEGER N */
/* DOUBLE PRECISION AMAX, SCOND */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION S( * ) */
/* COMPLEX*16 AP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQSP equilibrates a symmetric matrix A using the scaling factors */
/* > in the vector S. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > symmetric matrix A is stored. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the symmetric matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, the equilibrated matrix: diag(S) * A * diag(S), in */
/* > the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is DOUBLE PRECISION array, dimension (N) */
/* > The scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] SCOND */
/* > \verbatim */
/* > SCOND is DOUBLE PRECISION */
/* > Ratio of the smallest S(i) to the largest S(i). */
/* > \endverbatim */
/* > */
/* > \param[in] AMAX */
/* > \verbatim */
/* > AMAX is DOUBLE PRECISION */
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
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
void zlaqsp_(char *uplo, integer *n, doublecomplex *ap, doublereal *s, doublereal *scond,
             doublereal *amax, char *equed)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlaqsp inputs: uplo %c, n %" FLA_IS "", *uplo, *n);
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j, jc;
    doublereal cj, large;
    extern logical lsame_(char *, char *, integer, integer);
    doublereal small_val;
    extern doublereal dlamch_(char *);
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
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --s;
    --ap;
    /* Function Body */
    if(*n <= 0)
    {
        *(unsigned char *)equed = 'N';
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Initialize LARGE and SMALL. */
    small_val = dlamch_("Safe minimum") / dlamch_("Precision");
    large = 1. / small_val;
    if(*scond >= .1 && *amax >= small_val && *amax <= large)
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
            jc = 1;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                cj = s[j];
                i__2 = j;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = jc + i__ - 1;
                    d__1 = cj * s[i__];
                    i__4 = jc + i__ - 1;
                    z__1.r = d__1 * ap[i__4].r;
                    z__1.i = d__1 * ap[i__4].i; // , expr subst
                    ap[i__3].r = z__1.r;
                    ap[i__3].i = z__1.i; // , expr subst
                    /* L10: */
                }
                jc += j;
                /* L20: */
            }
        }
        else
        {
            /* Lower triangle of A is stored. */
            jc = 1;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                cj = s[j];
                i__2 = *n;
                for(i__ = j; i__ <= i__2; ++i__)
                {
                    i__3 = jc + i__ - j;
                    d__1 = cj * s[i__];
                    i__4 = jc + i__ - j;
                    z__1.r = d__1 * ap[i__4].r;
                    z__1.i = d__1 * ap[i__4].i; // , expr subst
                    ap[i__3].r = z__1.r;
                    ap[i__3].i = z__1.i; // , expr subst
                    /* L30: */
                }
                jc = jc + *n - j + 1;
                /* L40: */
            }
        }
        *(unsigned char *)equed = 'Y';
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLAQSP */
}
/* zlaqsp_ */
