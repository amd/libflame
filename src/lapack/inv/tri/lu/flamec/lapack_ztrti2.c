/* ztrti2.f -- translated by f2c (version 20190311).
   You must link the resulting object file with libf2c:
    on Microsoft Windows system, link with libf2c.lib;
    on Linux or Unix systems, link with .../path/to/libf2c.a -lm
    or, if you install libf2c.a in a standard place, with -lf2c -lm
    -- in that order, at the end of the command line, as in
        cc *.o -lf2c -lm
    Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

        http://www.netlib.org/f2c/libf2c.zip
*/
/*
 * Modifications Copyright (C) 2025 Advanced Micro Devices, Inc.  All
 * rights reserved.
 */
#include "FLA_f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTRTI2 computes the inverse of a triangular matrix (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTRTI2 + dependencies */
/* > <a
href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrti2.
f"> */
/* > [TGZ]</a> */
/* > <a
href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrti2.
f"> */
/* > [ZIP]</a> */
/* > <a
href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrti2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRTI2 computes the inverse of a complex upper or lower triangular */
/* > matrix. */
/* > */
/* > This is the Level 2 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the matrix A is upper or lower triangular. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          Specifies whether or not the matrix A is unit triangular. */
/* >          = 'N':  Non-unit triangular */
/* >          = 'U':  Unit triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the triangular matrix A.  If UPLO = 'U', the */
/* >          leading n by n upper triangular part of the array A contains */
/* >          the upper triangular matrix, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n by n lower triangular part of the array A contains */
/* >          the lower triangular matrix, and the strictly upper */
/* >          triangular part of A is not referenced.  If DIAG = 'U', the */
/* >          diagonal elements of A are also not referenced and are */
/* >          assumed to be 1. */
/* > */
/* >          On exit, the (triangular) inverse of the original matrix, in */
/* >          the same storage format. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup trti2 */


/* Z_DIV macro for complex division - follows z_div function logic */
#ifdef _WIN32
/* Windows implementation - manual complex division calculation */
#define Z_DIV_TRTI2(result, numerator, denominator)                                               \
    do                                                                                            \
    {                                                                                             \
        double temp;                                                                              \
        double s, xr_s, xi_s;                                                                     \
        s = fla_max(fabs((denominator)->real), fabs((denominator)->imag));                        \
        s = 1.0 / s;                                                                              \
        xr_s = (denominator)->real * s;                                                           \
        xi_s = (denominator)->imag * s;                                                           \
        temp = xr_s * (denominator)->real + xi_s * (denominator)->imag;                           \
        temp = 1.0 / temp;                                                                        \
        (result)->real = xr_s * temp;                                                             \
        (result)->imag = -xi_s * temp;                                                            \
        (result)->real = (result)->real * (numerator)->real + (result)->imag * (numerator)->imag; \
        (result)->imag = (result)->imag * (numerator)->real - (result)->real * (numerator)->imag; \
    } while(0)
#else
/* Non-Windows implementation - use C99 complex arithmetic */
#define Z_DIV_TRTI2(result, numerator, denominator)                                   \
    do                                                                                \
    {                                                                                 \
        dcomplex _a = *(numerator);                                                   \
        dcomplex _b = *(denominator);                                                 \
        double _Complex _ret_val = (_a.real + I * _a.imag) / (_b.real + I * _b.imag); \
        (result)->real = creal(_ret_val);                                             \
        (result)->imag = cimag(_ret_val);                                             \
    } while(0)
#endif

/*  ===================================================================== */
/* Subroutine */ void lapack_ztrti2(char *uplo, char *diag, integer *n, dcomplex *a, integer *lda,
                                    integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer j;
    dcomplex ajj;
    extern /* Subroutine */ int zscal_(integer *, dcomplex *, dcomplex *, integer *);
    extern int lsame_(char *, char *, integer a, integer b);
    logical upper;
    extern /* Subroutine */ int ztrmv_(char *, char *, char *, integer *, dcomplex *, integer *,
                                       dcomplex *, integer *),
        xerbla_(char *, integer *, ftnlen srname_len);
    logical nounit;

    /*  -- LAPACK computational routine -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(!nounit && !lsame_(diag, "U", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTRTI2", &i__1, (ftnlen)6);
        return;
    }

    if(upper)
    {

        /*        Compute inverse of upper triangular matrix. */

        i__1 = *n;
        if(nounit)
        {
            for(j = 1; j <= i__1; ++j)
            {
                dcomplex z__1;
                dcomplex c_b1 = {1.0, 0.0}; /* Complex constant 1+0i */
                i__2 = j + j * a_dim1;
                Z_DIV_TRTI2(&z__1, &c_b1, &a[j + j * a_dim1]);
                a[i__2].real = z__1.real, a[i__2].imag = z__1.imag;
                i__2 = j + j * a_dim1;
                z__1.real = -a[i__2].real, z__1.imag = -a[i__2].imag;
                ajj.real = z__1.real, ajj.imag = z__1.imag;

                /*           Compute elements 1:j-1 of j-th column. */

                i__2 = j - 1;
                ztrmv_("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &a[j * a_dim1 + 1],
                       &c__1);
                i__2 = j - 1;
                zscal_(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
            }
        }
        else
        {
            ajj.real = -1.;
            ajj.imag = 0.;
            for(j = 1; j <= i__1; ++j)
            {
                /*           Compute elements 1:j-1 of j-th column. */

                i__2 = j - 1;
                ztrmv_("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &a[j * a_dim1 + 1],
                       &c__1);
                i__2 = j - 1;
                zscal_(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
            }
        }
    }
    else
    {

        /*        Compute inverse of lower triangular matrix. */

        if(nounit)
        {
            for(j = *n; j >= 1; --j)
            {
                dcomplex z__1;
                dcomplex c_b1 = {1.0, 0.0}; /* Complex constant 1+0i */
                i__2 = j + j * a_dim1;
                Z_DIV_TRTI2(&z__1, &c_b1, &a[j + j * a_dim1]);
                a[i__2].real = z__1.real, a[i__2].imag = z__1.imag;
                i__2 = j + j * a_dim1;
                z__1.real = -a[i__2].real, z__1.imag = -a[i__2].imag;
                ajj.real = z__1.real, ajj.imag = z__1.imag;
                if(j < *n)
                {

                    /*              Compute elements j+1:n of j-th column. */

                    i__1 = *n - j;
                    ztrmv_("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 1) * a_dim1], lda,
                           &a[j + 1 + j * a_dim1], &c__1);
                    i__1 = *n - j;
                    zscal_(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
                }
            }
        }
        else
        {
            ajj.real = -1.;
            ajj.imag = 0.;
            for(j = *n; j >= 1; --j)
            {
                if(j < *n)
                {

                    /*              Compute elements j+1:n of j-th column. */

                    i__1 = *n - j;
                    ztrmv_("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 1) * a_dim1], lda,
                           &a[j + 1 + j * a_dim1], &c__1);
                    i__1 = *n - j;
                    zscal_(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
                }
            }
        }
    }

    return;

    /*     End of ZTRTI2 */

} /* ztrti2_ */
