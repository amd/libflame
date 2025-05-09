/* ../netlib/stpttf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b STPTTF copies a triangular matrix from the standard packed format (TP) to the rectangular full packed format (TF). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STPTTF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpttf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpttf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpttf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STPTTF( TRANSR, UPLO, N, AP, ARF, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANSR, UPLO */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AP( 0: * ), ARF( 0: * ) */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPTTF copies a triangular matrix A from standard packed format (TP) */
/* > to rectangular full packed format (TF). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANSR */
/* > \verbatim */
/* > TRANSR is CHARACTER*1 */
/* > = 'N': ARF in Normal format is wanted;
 */
/* > = 'T': ARF in Conjugate-transpose format is wanted. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': A is upper triangular;
 */
/* > = 'L': A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is REAL array, dimension ( N*(N+1)/2 ), */
/* > On entry, the upper or lower triangular matrix A, packed */
/* > columnwise in a linear array. The j-th column of A is stored */
/* > in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[out] ARF */
/* > \verbatim */
/* > ARF is REAL array, dimension ( N*(N+1)/2 ), */
/* > On exit, the upper or lower triangular matrix A stored in */
/* > RFP format. For a further discussion see Notes below. */
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
/* > \date September 2012 */
/* > \ingroup realOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > We first consider Rectangular Full Packed (RFP) Format when N is */
/* > even. We give an example where N = 6. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 05 00 */
/* > 11 12 13 14 15 10 11 */
/* > 22 23 24 25 20 21 22 */
/* > 33 34 35 30 31 32 33 */
/* > 44 45 40 41 42 43 44 */
/* > 55 50 51 52 53 54 55 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(4:6,0:2) consists of */
/* > the transpose of the first three columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:2,0:2) consists of */
/* > the transpose of the last three columns of AP lower. */
/* > This covers the case N even and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > 03 04 05 33 43 53 */
/* > 13 14 15 00 44 54 */
/* > 23 24 25 10 11 55 */
/* > 33 34 35 20 21 22 */
/* > 00 44 45 30 31 32 */
/* > 01 11 55 40 41 42 */
/* > 02 12 22 50 51 52 */
/* > */
/* > Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* > RFP A RFP A */
/* > */
/* > 03 13 23 33 00 01 02 33 00 10 20 30 40 50 */
/* > 04 14 24 34 44 11 12 43 44 11 21 31 41 51 */
/* > 05 15 25 35 45 55 22 53 54 55 22 32 42 52 */
/* > */
/* > */
/* > We then consider Rectangular Full Packed (RFP) Format when N is */
/* > odd. We give an example where N = 5. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 00 */
/* > 11 12 13 14 10 11 */
/* > 22 23 24 20 21 22 */
/* > 33 34 30 31 32 33 */
/* > 44 40 41 42 43 44 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(3:4,0:1) consists of */
/* > the transpose of the first two columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:1,1:2) consists of */
/* > the transpose of the last two columns of AP lower. */
/* > This covers the case N odd and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > 02 03 04 00 33 43 */
/* > 12 13 14 10 11 44 */
/* > 22 23 24 20 21 22 */
/* > 00 33 34 30 31 32 */
/* > 01 11 44 40 41 42 */
/* > */
/* > Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > RFP A RFP A */
/* > */
/* > 02 12 22 00 01 00 10 20 30 40 50 */
/* > 03 13 23 33 11 33 11 21 31 41 51 */
/* > 04 14 24 34 44 43 44 22 32 42 52 */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void stpttf_(char *transr, char *uplo, integer *n, real *ap, real *arf, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("stpttf inputs: transr %c ,uplo %c ,n %" FLA_IS "", *transr, *uplo, *n);
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Local variables */
    integer i__, j, k, n1, n2, ij, jp, js, lda, ijp;
    logical normaltransr;
    extern logical lsame_(char *, char *, integer, integer);
    logical lower;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical nisodd;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
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
    *info = 0;
    normaltransr = lsame_(transr, "N", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    if(!normaltransr && !lsame_(transr, "T", 1, 1))
    {
        *info = -1;
    }
    else if(!lower && !lsame_(uplo, "U", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STPTTF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 1)
    {
        if(normaltransr)
        {
            arf[0] = ap[0];
        }
        else
        {
            arf[0] = ap[0];
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Size of array ARF(0:NT-1) */
    /* Set N1 and N2 depending on LOWER */
    if(lower)
    {
        n2 = *n / 2;
        n1 = *n - n2;
    }
    else
    {
        n1 = *n / 2;
        n2 = *n - n1;
    }
    /* If N is odd, set NISODD = .TRUE. */
    /* If N is even, set K = N/2 and NISODD = .FALSE. */
    /* set lda of ARF^C;
    ARF^C is (0:(N+1)/2-1,0:N-noe) */
    /* where noe = 0 if n is even, noe = 1 if n is odd */
    if(*n % 2 == 0)
    {
        k = *n / 2;
        nisodd = FALSE_;
        lda = *n + 1;
    }
    else
    {
        nisodd = TRUE_;
        lda = *n;
    }
    /* ARF^C has lda rows and n+1-noe cols */
    if(!normaltransr)
    {
        lda = (*n + 1) / 2;
    }
    /* start execution: there are eight cases */
    if(nisodd)
    {
        /* N is odd */
        if(normaltransr)
        {
            /* N is odd and TRANSR = 'N' */
            if(lower)
            {
                /* N is odd, TRANSR = 'N', and UPLO = 'L' */
                ijp = 0;
                jp = 0;
                i__1 = n2;
                for(j = 0; j <= i__1; ++j)
                {
                    i__2 = *n - 1;
                    for(i__ = j; i__ <= i__2; ++i__)
                    {
                        ij = i__ + jp;
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    jp += lda;
                }
                i__1 = n2 - 1;
                for(i__ = 0; i__ <= i__1; ++i__)
                {
                    i__2 = n2;
                    for(j = i__ + 1; j <= i__2; ++j)
                    {
                        ij = i__ + j * lda;
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                }
            }
            else
            {
                /* N is odd, TRANSR = 'N', and UPLO = 'U' */
                ijp = 0;
                i__1 = n1 - 1;
                for(j = 0; j <= i__1; ++j)
                {
                    ij = n2 + j;
                    i__2 = j;
                    for(i__ = 0; i__ <= i__2; ++i__)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                        ij += lda;
                    }
                }
                js = 0;
                i__1 = *n - 1;
                for(j = n1; j <= i__1; ++j)
                {
                    ij = js;
                    i__2 = js + j;
                    for(ij = js; ij <= i__2; ++ij)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    js += lda;
                }
            }
        }
        else
        {
            /* N is odd and TRANSR = 'T' */
            if(lower)
            {
                /* N is odd, TRANSR = 'T', and UPLO = 'L' */
                ijp = 0;
                i__1 = n2;
                for(i__ = 0; i__ <= i__1; ++i__)
                {
                    i__2 = *n * lda - 1;
                    i__3 = lda;
                    for(ij = i__ * (lda + 1); i__3 < 0 ? ij >= i__2 : ij <= i__2; ij += i__3)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                }
                js = 1;
                i__1 = n2 - 1;
                for(j = 0; j <= i__1; ++j)
                {
                    i__3 = js + n2 - j - 1;
                    for(ij = js; ij <= i__3; ++ij)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    js = js + lda + 1;
                }
            }
            else
            {
                /* N is odd, TRANSR = 'T', and UPLO = 'U' */
                ijp = 0;
                js = n2 * lda;
                i__1 = n1 - 1;
                for(j = 0; j <= i__1; ++j)
                {
                    i__3 = js + j;
                    for(ij = js; ij <= i__3; ++ij)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    js += lda;
                }
                i__1 = n1;
                for(i__ = 0; i__ <= i__1; ++i__)
                {
                    i__3 = i__ + (n1 + i__) * lda;
                    i__2 = lda;
                    for(ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += i__2)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                }
            }
        }
    }
    else
    {
        /* N is even */
        if(normaltransr)
        {
            /* N is even and TRANSR = 'N' */
            if(lower)
            {
                /* N is even, TRANSR = 'N', and UPLO = 'L' */
                ijp = 0;
                jp = 0;
                i__1 = k - 1;
                for(j = 0; j <= i__1; ++j)
                {
                    i__2 = *n - 1;
                    for(i__ = j; i__ <= i__2; ++i__)
                    {
                        ij = i__ + 1 + jp;
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    jp += lda;
                }
                i__1 = k - 1;
                for(i__ = 0; i__ <= i__1; ++i__)
                {
                    i__2 = k - 1;
                    for(j = i__; j <= i__2; ++j)
                    {
                        ij = i__ + j * lda;
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                }
            }
            else
            {
                /* N is even, TRANSR = 'N', and UPLO = 'U' */
                ijp = 0;
                i__1 = k - 1;
                for(j = 0; j <= i__1; ++j)
                {
                    ij = k + 1 + j;
                    i__2 = j;
                    for(i__ = 0; i__ <= i__2; ++i__)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                        ij += lda;
                    }
                }
                js = 0;
                i__1 = *n - 1;
                for(j = k; j <= i__1; ++j)
                {
                    ij = js;
                    i__2 = js + j;
                    for(ij = js; ij <= i__2; ++ij)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    js += lda;
                }
            }
        }
        else
        {
            /* N is even and TRANSR = 'T' */
            if(lower)
            {
                /* N is even, TRANSR = 'T', and UPLO = 'L' */
                ijp = 0;
                i__1 = k - 1;
                for(i__ = 0; i__ <= i__1; ++i__)
                {
                    i__2 = (*n + 1) * lda - 1;
                    i__3 = lda;
                    for(ij = i__ + (i__ + 1) * lda; i__3 < 0 ? ij >= i__2 : ij <= i__2; ij += i__3)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                }
                js = 0;
                i__1 = k - 1;
                for(j = 0; j <= i__1; ++j)
                {
                    i__3 = js + k - j - 1;
                    for(ij = js; ij <= i__3; ++ij)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    js = js + lda + 1;
                }
            }
            else
            {
                /* N is even, TRANSR = 'T', and UPLO = 'U' */
                ijp = 0;
                js = (k + 1) * lda;
                i__1 = k - 1;
                for(j = 0; j <= i__1; ++j)
                {
                    i__3 = js + j;
                    for(ij = js; ij <= i__3; ++ij)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                    js += lda;
                }
                i__1 = k - 1;
                for(i__ = 0; i__ <= i__1; ++i__)
                {
                    i__3 = i__ + (k + i__) * lda;
                    i__2 = lda;
                    for(ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += i__2)
                    {
                        arf[ij] = ap[ijp];
                        ++ijp;
                    }
                }
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of STPTTF */
}
/* stpttf_ */
