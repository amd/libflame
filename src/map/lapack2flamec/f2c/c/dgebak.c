/* ../netlib/dgebak.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DGEBAK */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEBAK + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebak.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebak.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebak.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOB, SIDE */
/* INTEGER IHI, ILO, INFO, LDV, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION SCALE( * ), V( LDV, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEBAK forms the right or left eigenvectors of a real general matrix */
/* > by backward transformation on the computed eigenvectors of the */
/* > balanced matrix output by DGEBAL. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is CHARACTER*1 */
/* > Specifies the type of backward transformation required: */
/* > = 'N', do nothing, return immediately;
 */
/* > = 'P', do backward transformation for permutation only;
 */
/* > = 'S', do backward transformation for scaling only;
 */
/* > = 'B', do backward transformations for both permutation and */
/* > scaling. */
/* > JOB must be the same as the argument JOB supplied to DGEBAL. */
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
/* > The integers ILO and IHI determined by DGEBAL. */
/* > 1 <= ILO <= IHI <= N, if N > 0;
ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION array, dimension (N) */
/* > Details of the permutation and scaling factors, as returned */
/* > by DGEBAL. */
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
/* > transformed, as returned by DHSEIN or DTREVC. */
/* > On exit, V is overwritten by the transformed eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. LDV >= fla_max(1,N). */
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
/* > \date November 2011 */
/* > \ingroup doubleGEcomputational */
/* ===================================================================== */
/* Subroutine */
void dgebak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, doublereal *scale,
             integer *m, doublereal *v, integer *ldv, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgebak inputs: job %c, side %c, n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS
                      ", m %" FLA_IS ", ldv %" FLA_IS "",
                      *job, *side, *n, *ilo, *ihi, *m, *ldv);
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    /* Local variables */
    integer i__, k;
    doublereal s;
    integer ii;
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
    /* Decode and Test the input parameters */
    /* Parameter adjustments */
    --scale;
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
    else if(*ilo < 1 || *ilo > fla_max(1, *n))
    {
        *info = -4;
    }
    else if(*ihi < fla_min(*ilo, *n) || *ihi > *n)
    {
        *info = -5;
    }
    else if(*m < 0)
    {
        *info = -7;
    }
    else if(*ldv < fla_max(1, *n))
    {
        *info = -9;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGEBAK", &i__1, (ftnlen)6);
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
        if(rightv)
        {
            i__1 = *ihi;
            for(i__ = *ilo; i__ <= i__1; ++i__)
            {
                s = scale[i__];
                dscal_(m, &s, &v[i__ + v_dim1], ldv);
                /* L10: */
            }
        }
        if(leftv)
        {
            i__1 = *ihi;
            for(i__ = *ilo; i__ <= i__1; ++i__)
            {
                s = 1. / scale[i__];
                dscal_(m, &s, &v[i__ + v_dim1], ldv);
                /* L20: */
            }
        }
    }
    /* Backward permutation */
    /* For I = ILO-1 step -1 until 1, */
    /* IHI+1 step 1 until N do -- */
L30:
    if(lsame_(job, "P", 1, 1) || lsame_(job, "B", 1, 1))
    {
        if(rightv)
        {
            i__1 = *n;
            for(ii = 1; ii <= i__1; ++ii)
            {
                i__ = ii;
                if(i__ >= *ilo && i__ <= *ihi)
                {
                    goto L40;
                }
                if(i__ < *ilo)
                {
                    i__ = *ilo - ii;
                }
                k = (integer)scale[i__];
                if(k == i__)
                {
                    goto L40;
                }
                dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
            L40:;
            }
        }
        if(leftv)
        {
            i__1 = *n;
            for(ii = 1; ii <= i__1; ++ii)
            {
                i__ = ii;
                if(i__ >= *ilo && i__ <= *ihi)
                {
                    goto L50;
                }
                if(i__ < *ilo)
                {
                    i__ = *ilo - ii;
                }
                k = (integer)scale[i__];
                if(k == i__)
                {
                    goto L50;
                }
                dswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
            L50:;
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DGEBAK */
}
/* dgebak_ */
