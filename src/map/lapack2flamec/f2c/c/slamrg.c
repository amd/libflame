/* ../netlib/slamrg.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLAMRG creates a permutation list to merge the entries of two independently sorted sets into a single set sorted in ascending order. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAMRG + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slamrg.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slamrg.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slamrg.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAMRG( N1, N2, A, STRD1, STRD2, INDEX ) */
/* .. Scalar Arguments .. */
/* INTEGER N1, N2, STRD1, STRD2 */
/* .. */
/* .. Array Arguments .. */
/* INTEGER INDEX( * ) */
/* REAL A( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAMRG will create a permutation list which will merge the elements */
/* > of A (which is composed of two independently sorted sets) into a */
/* > single set which is sorted in ascending order. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N1 */
/* > \verbatim */
/* > N1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* > N2 is INTEGER */
/* > These arguements contain the respective lengths of the two */
/* > sorted lists to be merged. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (N1+N2) */
/* > The first N1 elements of A contain a list of numbers which */
/* > are sorted in either ascending or descending order. Likewise */
/* > for the final N2 elements. */
/* > \endverbatim */
/* > */
/* > \param[in] STRD1 */
/* > \verbatim */
/* > STRD1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] STRD2 */
/* > \verbatim */
/* > STRD2 is INTEGER */
/* > These are the strides to be taken through the array A. */
/* > Allowable strides are 1 and -1. They indicate whether a */
/* > subset of A is sorted in ascending (STRDx = 1) or descending */
/* > (STRDx = -1) order. */
/* > \endverbatim */
/* > */
/* > \param[out] INDEX */
/* > \verbatim */
/* > INDEX is INTEGER array, dimension (N1+N2) */
/* > On exit this array will contain a permutation such that */
/* > if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be */
/* > sorted in ascending order. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void slamrg_(aocl_int_t *n1, aocl_int_t *n2, real *a, aocl_int_t *strd1, aocl_int_t *strd2,
             aocl_int_t *index)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_slamrg(n1, n2, a, strd1, strd2, index);
#else
    aocl_int64_t n1_64 = *n1;
    aocl_int64_t n2_64 = *n2;
    aocl_int64_t strd1_64 = *strd1;
    aocl_int64_t strd2_64 = *strd2;

    aocl_lapack_slamrg(&n1_64, &n2_64, a, &strd1_64, &strd2_64, index);
#endif
}

void aocl_lapack_slamrg(aocl_int64_t *n1, aocl_int64_t *n2, real *a, aocl_int64_t *strd1,
                        aocl_int64_t *strd2, aocl_int_t *index)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slamrg inputs: n1 %" FLA_IS ", n2 %" FLA_IS ", strd1 %" FLA_IS
                      ", strd2 %" FLA_IS "",
                      *n1, *n2, *strd1, *strd2);
    /* System generated locals */
    aocl_int64_t i__1;
    /* Local variables */
    aocl_int64_t i__, ind1, ind2, n1sv, n2sv;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --index;
    --a;
    /* Function Body */
    n1sv = *n1;
    n2sv = *n2;
    if(*strd1 > 0)
    {
        ind1 = 1;
    }
    else
    {
        ind1 = *n1;
    }
    if(*strd2 > 0)
    {
        ind2 = *n1 + 1;
    }
    else
    {
        ind2 = *n1 + *n2;
    }
    i__ = 1;
    /* while ( (N1SV > 0) & (N2SV > 0) ) */
L10:
    if(n1sv > 0 && n2sv > 0)
    {
        if(a[ind1] <= a[ind2])
        {
            index[i__] = (aocl_int_t)(ind1);
            ++i__;
            ind1 += *strd1;
            --n1sv;
        }
        else
        {
            index[i__] = (aocl_int_t)(ind2);
            ++i__;
            ind2 += *strd2;
            --n2sv;
        }
        goto L10;
    }
    /* end while */
    if(n1sv == 0)
    {
        i__1 = n2sv;
        for(n1sv = 1; n1sv <= i__1; ++n1sv)
        {
            index[i__] = (aocl_int_t)(ind2);
            ++i__;
            ind2 += *strd2;
            /* L20: */
        }
    }
    else
    {
        /* N2SV .EQ. 0 */
        i__1 = n1sv;
        for(n2sv = 1; n2sv <= i__1; ++n2sv)
        {
            index[i__] = (aocl_int_t)(ind1);
            ++i__;
            ind1 += *strd1;
            /* L30: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLAMRG */
}
/* slamrg_ */
