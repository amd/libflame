/* ../netlib/dggbak.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DGGBAK */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGGBAK + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbak.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbak.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbak.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, */
/* LDV, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOB, SIDE */
/* INTEGER IHI, ILO, INFO, LDV, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION LSCALE( * ), RSCALE( * ), V( LDV, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGBAK forms the right or left eigenvectors of a real generalized */
/* > eigenvalue problem A*x = lambda*B*x, by backward transformation on */
/* > the computed eigenvectors of the balanced pair of matrices output by */
/* > DGGBAL. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is CHARACTER*1 */
/* > Specifies the type of backward transformation required: */
/* > = 'N': do nothing, return immediately;
 */
/* > = 'P': do backward transformation for permutation only;
 */
/* > = 'S': do backward transformation for scaling only;
 */
/* > = 'B': do backward transformations for both permutation and */
/* > scaling. */
/* > JOB must be the same as the argument JOB supplied to DGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'R': V contains right eigenvectors;
 */
/* > = 'L': V contains left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of rows of the matrix V. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > The integers ILO and IHI determined by DGGBAL. */
/* > 1 <= ILO <= IHI <= N, if N > 0;
ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] LSCALE */
/* > \verbatim */
/* > LSCALE is DOUBLE PRECISION array, dimension (N) */
/* > Details of the permutations and/or scaling factors applied */
/* > to the left side of A and B, as returned by DGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] RSCALE */
/* > \verbatim */
/* > RSCALE is DOUBLE PRECISION array, dimension (N) */
/* > Details of the permutations and/or scaling factors applied */
/* > to the right side of A and B, as returned by DGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of columns of the matrix V. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension (LDV,M) */
/* > On entry, the matrix of right or left eigenvectors to be */
/* > transformed, as returned by DTGEVC. */
/* > On exit, V is overwritten by the transformed eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the matrix V. LDV >= fla_max(1,N). */
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
/* > \date December 2016 */
/* > \ingroup doubleGBcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > See R.C. Ward, Balancing the generalized eigenvalue problem, */
/* > SIAM J. Sci. Stat. Comp. 2 (1981), 141-152. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dggbak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, doublereal *lscale,
             doublereal *rscale, integer *m, doublereal *v, integer *ldv, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dggbak inputs: job %c, side %c, n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS
                      ", m %" FLA_IS ", ldv %" FLA_IS "",
                      *job, *side, *n, *ilo, *ihi, *m, *ldv);
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    /* Local variables */
    integer i__, k;
    extern /* Subroutine */
        void
        dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical leftv;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical rightv;
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    --lscale;
    --rscale;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    /* Function Body */
    rightv = lsame_(side, "R", 1, 1);
    leftv = lsame_(side, "L", 1, 1);
    *info = 0;
    if(!lsame_(job, "N", 1, 1) && !lsame_(job, "P", 1, 1) && !lsame_(job, "S", 1, 1)
       && !lsame_(job, "B", 1, 1))
    {
        *info = -1;
    }
    else if(!rightv && !leftv)
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*ilo < 1)
    {
        *info = -4;
    }
    else if(*n == 0 && *ihi == 0 && *ilo != 1)
    {
        *info = -4;
    }
    else if(*n > 0 && (*ihi < *ilo || *ihi > fla_max(1, *n)))
    {
        *info = -5;
    }
    else if(*n == 0 && *ilo == 1 && *ihi != 0)
    {
        *info = -5;
    }
    else if(*m < 0)
    {
        *info = -8;
    }
    else if(*ldv < fla_max(1, *n))
    {
        *info = -10;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGGBAK", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*m == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(lsame_(job, "N", 1, 1))
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*ilo == *ihi)
    {
        goto L30;
    }
    /* Backward balance */
    if(lsame_(job, "S", 1, 1) || lsame_(job, "B", 1, 1))
    {
        /* Backward transformation on right eigenvectors */
        if(rightv)
        {
            i__1 = *ihi;
            for(i__ = *ilo; i__ <= i__1; ++i__)
            {
                dscal_(m, &rscale[i__], &v[i__ + v_dim1], ldv);
                /* L10: */
            }
        }
        /* Backward transformation on left eigenvectors */
        if(leftv)
        {
            i__1 = *ihi;
            for(i__ = *ilo; i__ <= i__1; ++i__)
            {
                dscal_(m, &lscale[i__], &v[i__ + v_dim1], ldv);
                /* L20: */
            }
        }
    }
    /* Backward permutation */
L30:
    if(lsame_(job, "P", 1, 1) || lsame_(job, "B", 1, 1))
    {
        /* Backward permutation on right eigenvectors */
        if(rightv)
        {
            if(*ilo == 1)
            {
                goto L50;
            }
            for(i__ = *ilo - 1; i__ >= 1; --i__)
            {
                k = (integer)rscale[i__];
                if(k == i__)
                {
                    goto L40;
                }
                dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
            L40:;
            }
        L50:
            if(*ihi == *n)
            {
                goto L70;
            }
            i__1 = *n;
            for(i__ = *ihi + 1; i__ <= i__1; ++i__)
            {
                k = (integer)rscale[i__];
                if(k == i__)
                {
                    goto L60;
                }
                dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
            L60:;
            }
        }
        /* Backward permutation on left eigenvectors */
    L70:
        if(leftv)
        {
            if(*ilo == 1)
            {
                goto L90;
            }
            for(i__ = *ilo - 1; i__ >= 1; --i__)
            {
                k = (integer)lscale[i__];
                if(k == i__)
                {
                    goto L80;
                }
                dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
            L80:;
            }
        L90:
            if(*ihi == *n)
            {
                goto L110;
            }
            i__1 = *n;
            for(i__ = *ihi + 1; i__ <= i__1; ++i__)
            {
                k = (integer)lscale[i__];
                if(k == i__)
                {
                    goto L100;
                }
                dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
            L100:;
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
L110:
    return;
    /* End of DGGBAK */
}
/* dggbak_ */
