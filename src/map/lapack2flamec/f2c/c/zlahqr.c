/* ./zlahqr.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c__2 = 2;
/* > \brief \b ZLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg
 * matrix, using th e double-shift/single-shift QR algorithm. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAHQR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahqr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahqr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahqr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, */
/* IHIZ, Z, LDZ, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/* LOGICAL WANTT, WANTZ */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 H( LDH, * ), W( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAHQR is an auxiliary routine called by CHSEQR to update the */
/* > eigenvalues and Schur decomposition already computed by CHSEQR, by */
/* > dealing with the Hessenberg submatrix in rows and columns ILO to */
/* > IHI. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTT */
/* > \verbatim */
/* > WANTT is LOGICAL */
/* > = .TRUE. : the full Schur form T is required;
 */
/* > = .FALSE.: only eigenvalues are required. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > = .TRUE. : the matrix of Schur vectors Z is required;
 */
/* > = .FALSE.: Schur vectors are not required. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix H. N >= 0. */
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
/* > It is assumed that H is already upper triangular in rows and */
/* > columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1). */
/* > ZLAHQR works primarily with the Hessenberg submatrix in rows */
/* > and columns ILO to IHI, but applies transformations to all of */
/* > H if WANTT is .TRUE.. */
/* > 1 <= ILO <= fla_max(1,IHI);
IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* > H is COMPLEX*16 array, dimension (LDH,N) */
/* > On entry, the upper Hessenberg matrix H. */
/* > On exit, if INFO is zero and if WANTT is .TRUE., then H */
/* > is upper triangular in rows and columns ILO:IHI. If INFO */
/* > is zero and if WANTT is .FALSE., then the contents of H */
/* > are unspecified on exit. The output state of H in case */
/* > INF is positive is below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > The leading dimension of the array H. LDH >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX*16 array, dimension (N) */
/* > The computed eigenvalues ILO to IHI are stored in the */
/* > corresponding elements of W. If WANTT is .TRUE., the */
/* > eigenvalues are stored in the same order as on the diagonal */
/* > of the Schur form returned in H, with W(i) = H(i,i). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* > ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* > IHIZ is INTEGER */
/* > Specify the rows of Z to which transformations must be */
/* > applied if WANTZ is .TRUE.. */
/* > 1 <= ILOZ <= ILO;
IHI <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX*16 array, dimension (LDZ,N) */
/* > If WANTZ is .TRUE., on entry Z must contain the current */
/* > matrix Z of transformations accumulated by CHSEQR, and on */
/* > exit Z has been updated;
transformations are applied only to */
/* > the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
/* > If WANTZ is .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > > 0: if INFO = i, ZLAHQR failed to compute all the */
/* > eigenvalues ILO to IHI in a total of 30 iterations */
/* > per eigenvalue;
elements i+1:ihi of W contain */
/* > those eigenvalues which have been successfully */
/* > computed. */
/* > */
/* > If INFO > 0 and WANTT is .FALSE., then on exit, */
/* > the remaining unconverged eigenvalues are the */
/* > eigenvalues of the upper Hessenberg matrix */
/* > rows and columns ILO through INFO of the final, */
/* > output value of H. */
/* > */
/* > If INFO > 0 and WANTT is .TRUE., then on exit */
/* > (*) (initial value of H)*U = U*(final value of H) */
/* > where U is an orthogonal matrix. The final */
/* > value of H is upper Hessenberg and triangular in */
/* > rows and columns INFO+1 through IHI. */
/* > */
/* > If INFO > 0 and WANTZ is .TRUE., then on exit */
/* > (final value of Z) = (initial value of Z)*U */
/* > where U is the orthogonal matrix in (*) */
/* > (regardless of the value of WANTT.) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup lahqr */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > 02-96 Based on modifications by */
/* > David Day, Sandia National Laboratory, USA */
/* > */
/* > 12-04 Further modifications by */
/* > Ralph Byers, University of Kansas, USA */
/* > This is a modified version of ZLAHQR from LAPACK version 3.0. */
/* > It is (1) more robust against overflow and underflow and */
/* > (2) adopts the more conservative Ahues & Tisseur stopping */
/* > criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlahqr_(logical *wantt, logical *wantz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi,
             dcomplex *h__, aocl_int_t *ldh, dcomplex *w, aocl_int_t *iloz,
             aocl_int_t *ihiz, dcomplex *z__, aocl_int_t *ldz, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlahqr(wantt, wantz, n, ilo, ihi, h__, ldh, w, iloz, ihiz, z__, ldz, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ilo_64 = *ilo;
    aocl_int64_t ihi_64 = *ihi;
    aocl_int64_t ldh_64 = *ldh;
    aocl_int64_t iloz_64 = *iloz;
    aocl_int64_t ihiz_64 = *ihiz;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zlahqr(wantt, wantz, &n_64, &ilo_64, &ihi_64, h__, &ldh_64, w, &iloz_64, &ihiz_64,
                       z__, &ldz_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zlahqr(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo,
                        aocl_int64_t *ihi, dcomplex *h__, aocl_int64_t *ldh, dcomplex *w,
                        aocl_int64_t *iloz, aocl_int64_t *ihiz, dcomplex *z__,
                        aocl_int64_t *ldz, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlahqr inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", ldh %" FLA_IS
                      ", iloz %" FLA_IS ", ihiz %" FLA_IS ", ldz %" FLA_IS "",
                      *n, *ilo, *ihi, *ldh, *iloz, *ihiz, *ldz);
    /* System generated locals */
    aocl_int64_t h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    dcomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;
    /* Builtin functions */
    double d_imag(dcomplex *);
    void d_cnjg(dcomplex *, dcomplex *);
    double z_abs(dcomplex *);
    void z_sqrt(dcomplex *, dcomplex *),
        pow_zi(dcomplex *, dcomplex *, aocl_int64_t *);
    /* Local variables */
    aocl_int64_t i__, j, k, l, m;
    doublereal s;
    dcomplex t, u, v[2], x, y;
    aocl_int64_t i1, i2;
    dcomplex t1;
    doublereal t2;
    dcomplex v2;
    doublereal aa, ab, ba, bb, h10;
    dcomplex h11;
    doublereal h21;
    dcomplex h22, sc;
    aocl_int64_t nh, nz;
    doublereal sx;
    aocl_int64_t jhi;
    dcomplex h11s;
    aocl_int64_t jlo, its;
    doublereal ulp;
    dcomplex sum;
    doublereal tst;
    dcomplex temp;
    aocl_int64_t kdefl;
    aocl_int64_t itmax;
    doublereal rtemp;
    extern doublereal dlamch_(char *);
    doublereal safmin;
    extern /* Double Complex */
        void
        zladiv_f2c_(dcomplex *, dcomplex *, dcomplex *);
    doublereal smlnum;
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ========================================================= */
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    /* Function Body */
    *info = 0;
    i2 = 0;
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*ilo == *ihi)
    {
        i__1 = *ilo;
        i__2 = *ilo + *ilo * h_dim1;
        w[i__1].real = h__[i__2].real;
        w[i__1].imag = h__[i__2].imag; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ==== clear out the trash ==== */
    i__1 = *ihi - 3;
    for(j = *ilo; j <= i__1; ++j)
    {
        i__2 = j + 2 + j * h_dim1;
        h__[i__2].real = 0.;
        h__[i__2].imag = 0.; // , expr subst
        i__2 = j + 3 + j * h_dim1;
        h__[i__2].real = 0.;
        h__[i__2].imag = 0.; // , expr subst
        /* L10: */
    }
    if(*ilo <= *ihi - 2)
    {
        i__1 = *ihi + (*ihi - 2) * h_dim1;
        h__[i__1].real = 0.;
        h__[i__1].imag = 0.; // , expr subst
    }
    /* ==== ensure that subdiagonal entries are real ==== */
    if(*wantt)
    {
        jlo = 1;
        jhi = *n;
    }
    else
    {
        jlo = *ilo;
        jhi = *ihi;
    }
    i__1 = *ihi;
    for(i__ = *ilo + 1; i__ <= i__1; ++i__)
    {
        if(d_imag(&h__[i__ + (i__ - 1) * h_dim1]) != 0.)
        {
            /* ==== The following redundant normalization */
            /* . avoids problems with both gradual and */
            /* . sudden underflow in ABS(H(I,I-1)) ==== */
            i__2 = i__ + (i__ - 1) * h_dim1;
            i__3 = i__ + (i__ - 1) * h_dim1;
            d__3 = (d__1 = h__[i__3].real, f2c_dabs(d__1))
                   + (d__2 = d_imag(&h__[i__ + (i__ - 1) * h_dim1]), f2c_dabs(d__2));
            z__1.real = h__[i__2].real / d__3;
            z__1.imag = h__[i__2].imag / d__3; // , expr subst
            sc.real = z__1.real;
            sc.imag = z__1.imag; // , expr subst
            d_cnjg(&z__2, &sc);
            d__1 = z_abs(&sc);
            z__1.real = z__2.real / d__1;
            z__1.imag = z__2.imag / d__1; // , expr subst
            sc.real = z__1.real;
            sc.imag = z__1.imag; // , expr subst
            i__2 = i__ + (i__ - 1) * h_dim1;
            d__1 = z_abs(&h__[i__ + (i__ - 1) * h_dim1]);
            h__[i__2].real = d__1;
            h__[i__2].imag = 0.; // , expr subst
            i__2 = jhi - i__ + 1;
            aocl_blas_zscal(&i__2, &sc, &h__[i__ + i__ * h_dim1], ldh);
            /* Computing MIN */
            i__3 = jhi;
            i__4 = i__ + 1; // , expr subst
            i__2 = fla_min(i__3, i__4) - jlo + 1;
            d_cnjg(&z__1, &sc);
            aocl_blas_zscal(&i__2, &z__1, &h__[jlo + i__ * h_dim1], &c__1);
            if(*wantz)
            {
                i__2 = *ihiz - *iloz + 1;
                d_cnjg(&z__1, &sc);
                aocl_blas_zscal(&i__2, &z__1, &z__[*iloz + i__ * z_dim1], &c__1);
            }
        }
        /* L20: */
    }
    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;
    /* Set machine-dependent constants for the stopping criterion. */
    safmin = dlamch_("SAFE MINIMUM");
    ulp = dlamch_("PRECISION");
    smlnum = safmin * ((doublereal)nh / ulp);
    /* I1 and I2 are the indices of the first row and last column of H */
    /* to which transformations must be applied. If eigenvalues only are */
    /* being computed, I1 and I2 are set inside the main loop. */
    if(*wantt)
    {
        i1 = 1;
        i2 = *n;
    }
    /* ITMAX is the total number of QR iterations allowed. */
    itmax = fla_max(10, nh) * 30;
    /* KDEFL counts the number of iterations since a deflation */
    kdefl = 0;
    /* The main loop begins here. I is the loop index and decreases from */
    /* IHI to ILO in steps of 1. Each iteration of the loop works */
    /* with the active submatrix in rows and columns L to I. */
    /* Eigenvalues I+1 to IHI have already converged. Either L = ILO, or */
    /* H(L,L-1) is negligible so that the matrix splits. */
    i__ = *ihi;
L30:
    if(i__ < *ilo)
    {
        goto L150;
    }
    /* Perform QR iterations on rows and columns ILO to I until a */
    /* submatrix of order 1 splits off at the bottom because a */
    /* subdiagonal element has become negligible. */
    l = *ilo;
    i__1 = itmax;
    for(its = 0; its <= i__1; ++its)
    {
        /* Look for a single small subdiagonal element. */
        i__2 = l + 1;
        for(k = i__; k >= i__2; --k)
        {
            i__3 = k + (k - 1) * h_dim1;
            if((d__1 = h__[i__3].real, f2c_dabs(d__1))
                   + (d__2 = d_imag(&h__[k + (k - 1) * h_dim1]), f2c_dabs(d__2))
               <= smlnum)
            {
                goto L50;
            }
            i__3 = k - 1 + (k - 1) * h_dim1;
            i__4 = k + k * h_dim1;
            tst = (d__1 = h__[i__3].real, f2c_dabs(d__1))
                  + (d__2 = d_imag(&h__[k - 1 + (k - 1) * h_dim1]), f2c_dabs(d__2))
                  + ((d__3 = h__[i__4].real, f2c_dabs(d__3))
                     + (d__4 = d_imag(&h__[k + k * h_dim1]), f2c_dabs(d__4)));
            if(tst == 0.)
            {
                if(k - 2 >= *ilo)
                {
                    i__3 = k - 1 + (k - 2) * h_dim1;
                    tst += (d__1 = h__[i__3].real, f2c_dabs(d__1));
                }
                if(k + 1 <= *ihi)
                {
                    i__3 = k + 1 + k * h_dim1;
                    tst += (d__1 = h__[i__3].real, f2c_dabs(d__1));
                }
            }
            /* ==== The following is a conservative small subdiagonal */
            /* . deflation criterion due to Ahues & Tisseur (LAWN 122, */
            /* . 1997). It has better mathematical foundation and */
            /* . improves accuracy in some examples. ==== */
            i__3 = k + (k - 1) * h_dim1;
            if((d__1 = h__[i__3].real, f2c_dabs(d__1)) <= ulp * tst)
            {
                /* Computing MAX */
                i__3 = k + (k - 1) * h_dim1;
                i__4 = k - 1 + k * h_dim1;
                d__5 = (d__1 = h__[i__3].real, f2c_dabs(d__1))
                       + (d__2 = d_imag(&h__[k + (k - 1) * h_dim1]), f2c_dabs(d__2));
                d__6 = (d__3 = h__[i__4].real, f2c_dabs(d__3))
                       + (d__4 = d_imag(&h__[k - 1 + k * h_dim1]), f2c_dabs(d__4)); // , expr subst
                ab = fla_max(d__5, d__6);
                /* Computing MIN */
                i__3 = k + (k - 1) * h_dim1;
                i__4 = k - 1 + k * h_dim1;
                d__5 = (d__1 = h__[i__3].real, f2c_dabs(d__1))
                       + (d__2 = d_imag(&h__[k + (k - 1) * h_dim1]), f2c_dabs(d__2));
                d__6 = (d__3 = h__[i__4].real, f2c_dabs(d__3))
                       + (d__4 = d_imag(&h__[k - 1 + k * h_dim1]), f2c_dabs(d__4)); // , expr subst
                ba = fla_min(d__5, d__6);
                i__3 = k - 1 + (k - 1) * h_dim1;
                i__4 = k + k * h_dim1;
                z__2.real = h__[i__3].real - h__[i__4].real;
                z__2.imag = h__[i__3].imag - h__[i__4].imag; // , expr subst
                z__1.real = z__2.real;
                z__1.imag = z__2.imag; // , expr subst
                /* Computing MAX */
                i__5 = k + k * h_dim1;
                d__5 = (d__1 = h__[i__5].real, f2c_dabs(d__1))
                       + (d__2 = d_imag(&h__[k + k * h_dim1]), f2c_dabs(d__2));
                d__6 = (d__3 = z__1.real, f2c_dabs(d__3))
                       + (d__4 = d_imag(&z__1), f2c_dabs(d__4)); // , expr subst
                aa = fla_max(d__5, d__6);
                i__3 = k - 1 + (k - 1) * h_dim1;
                i__4 = k + k * h_dim1;
                z__2.real = h__[i__3].real - h__[i__4].real;
                z__2.imag = h__[i__3].imag - h__[i__4].imag; // , expr subst
                z__1.real = z__2.real;
                z__1.imag = z__2.imag; // , expr subst
                /* Computing MIN */
                i__5 = k + k * h_dim1;
                d__5 = (d__1 = h__[i__5].real, f2c_dabs(d__1))
                       + (d__2 = d_imag(&h__[k + k * h_dim1]), f2c_dabs(d__2));
                d__6 = (d__3 = z__1.real, f2c_dabs(d__3))
                       + (d__4 = d_imag(&z__1), f2c_dabs(d__4)); // , expr subst
                bb = fla_min(d__5, d__6);
                s = aa + ab;
                /* Computing MAX */
                d__1 = smlnum;
                d__2 = ulp * (bb * (aa / s)); // , expr subst
                if(ba * (ab / s) <= fla_max(d__1, d__2))
                {
                    goto L50;
                }
            }
            /* L40: */
        }
    L50:
        l = k;
        if(l > *ilo)
        {
            /* H(L,L-1) is negligible */
            i__2 = l + (l - 1) * h_dim1;
            h__[i__2].real = 0.;
            h__[i__2].imag = 0.; // , expr subst
        }
        /* Exit from loop if a submatrix of order 1 has split off. */
        if(l >= i__)
        {
            goto L140;
        }
        ++kdefl;
        /* Now the active submatrix is in rows and columns L to I. If */
        /* eigenvalues only are being computed, only the active submatrix */
        /* need be transformed. */
        if(!(*wantt))
        {
            i1 = l;
            i2 = i__;
        }
        if(kdefl % 20 == 0)
        {
            /* Exceptional shift. */
            i__2 = i__ + (i__ - 1) * h_dim1;
            s = (d__1 = h__[i__2].real, f2c_dabs(d__1)) * .75;
            i__2 = i__ + i__ * h_dim1;
            z__1.real = s + h__[i__2].real;
            z__1.imag = h__[i__2].imag; // , expr subst
            t.real = z__1.real;
            t.imag = z__1.imag; // , expr subst
        }
        else if(kdefl % 10 == 0)
        {
            /* Exceptional shift. */
            i__2 = l + 1 + l * h_dim1;
            s = (d__1 = h__[i__2].real, f2c_dabs(d__1)) * .75;
            i__2 = l + l * h_dim1;
            z__1.real = s + h__[i__2].real;
            z__1.imag = h__[i__2].imag; // , expr subst
            t.real = z__1.real;
            t.imag = z__1.imag; // , expr subst
        }
        else
        {
            /* Wilkinson's shift. */
            i__2 = i__ + i__ * h_dim1;
            t.real = h__[i__2].real;
            t.imag = h__[i__2].imag; // , expr subst
            z_sqrt(&z__2, &h__[i__ - 1 + i__ * h_dim1]);
            z_sqrt(&z__3, &h__[i__ + (i__ - 1) * h_dim1]);
            z__1.real = z__2.real * z__3.real - z__2.imag * z__3.imag;
            z__1.imag = z__2.real * z__3.imag + z__2.imag * z__3.real; // , expr subst
            u.real = z__1.real;
            u.imag = z__1.imag; // , expr subst
            s = (d__1 = u.real, f2c_dabs(d__1)) + (d__2 = d_imag(&u), f2c_dabs(d__2));
            if(s != 0.)
            {
                i__2 = i__ - 1 + (i__ - 1) * h_dim1;
                z__2.real = h__[i__2].real - t.real;
                z__2.imag = h__[i__2].imag - t.imag; // , expr subst
                z__1.real = z__2.real * .5;
                z__1.imag = z__2.imag * .5; // , expr subst
                x.real = z__1.real;
                x.imag = z__1.imag; // , expr subst
                sx = (d__1 = x.real, f2c_dabs(d__1)) + (d__2 = d_imag(&x), f2c_dabs(d__2));
                /* Computing MAX */
                d__3 = s;
                d__4 = (d__1 = x.real, f2c_dabs(d__1))
                       + (d__2 = d_imag(&x), f2c_dabs(d__2)); // , expr subst
                s = fla_max(d__3, d__4);
                z__5.real = x.real / s;
                z__5.imag = x.imag / s; // , expr subst
                pow_zi(&z__4, &z__5, &c__2);
                z__7.real = u.real / s;
                z__7.imag = u.imag / s; // , expr subst
                pow_zi(&z__6, &z__7, &c__2);
                z__3.real = z__4.real + z__6.real;
                z__3.imag = z__4.imag + z__6.imag; // , expr subst
                z_sqrt(&z__2, &z__3);
                z__1.real = s * z__2.real;
                z__1.imag = s * z__2.imag; // , expr subst
                y.real = z__1.real;
                y.imag = z__1.imag; // , expr subst
                if(sx > 0.)
                {
                    z__1.real = x.real / sx;
                    z__1.imag = x.imag / sx; // , expr subst
                    z__2.real = x.real / sx;
                    z__2.imag = x.imag / sx; // , expr subst
                    if(z__1.real * y.real + d_imag(&z__2) * d_imag(&y) < 0.)
                    {
                        z__3.real = -y.real;
                        z__3.imag = -y.imag; // , expr subst
                        y.real = z__3.real;
                        y.imag = z__3.imag; // , expr subst
                    }
                }
                z__4.real = x.real + y.real;
                z__4.imag = x.imag + y.imag; // , expr subst
                zladiv_f2c_(&z__3, &u, &z__4);
                z__2.real = u.real * z__3.real - u.imag * z__3.imag;
                z__2.imag = u.real * z__3.imag + u.imag * z__3.real; // , expr subst
                z__1.real = t.real - z__2.real;
                z__1.imag = t.imag - z__2.imag; // , expr subst
                t.real = z__1.real;
                t.imag = z__1.imag; // , expr subst
            }
        }
        /* Look for two consecutive small subdiagonal elements. */
        i__2 = l + 1;
        for(m = i__ - 1; m >= i__2; --m)
        {
            /* Determine the effect of starting the single-shift QR */
            /* iteration at row M, and see if this would make H(M,M-1) */
            /* negligible. */
            i__3 = m + m * h_dim1;
            h11.real = h__[i__3].real;
            h11.imag = h__[i__3].imag; // , expr subst
            i__3 = m + 1 + (m + 1) * h_dim1;
            h22.real = h__[i__3].real;
            h22.imag = h__[i__3].imag; // , expr subst
            z__1.real = h11.real - t.real;
            z__1.imag = h11.imag - t.imag; // , expr subst
            h11s.real = z__1.real;
            h11s.imag = z__1.imag; // , expr subst
            i__3 = m + 1 + m * h_dim1;
            h21 = h__[i__3].real;
            s = (d__1 = h11s.real, f2c_dabs(d__1)) + (d__2 = d_imag(&h11s), f2c_dabs(d__2))
                + f2c_dabs(h21);
            z__1.real = h11s.real / s;
            z__1.imag = h11s.imag / s; // , expr subst
            h11s.real = z__1.real;
            h11s.imag = z__1.imag; // , expr subst
            h21 /= s;
            v[0].real = h11s.real;
            v[0].imag = h11s.imag; // , expr subst
            v[1].real = h21;
            v[1].imag = 0.; // , expr subst
            i__3 = m + (m - 1) * h_dim1;
            h10 = h__[i__3].real;
            if(f2c_dabs(h10) * f2c_dabs(h21)
               <= ulp
                      * (((d__1 = h11s.real, f2c_dabs(d__1)) + (d__2 = d_imag(&h11s), f2c_dabs(d__2)))
                         * ((d__3 = h11.real, f2c_dabs(d__3)) + (d__4 = d_imag(&h11), f2c_dabs(d__4))
                            + ((d__5 = h22.real, f2c_dabs(d__5))
                               + (d__6 = d_imag(&h22), f2c_dabs(d__6))))))
            {
                goto L70;
            }
            /* L60: */
        }
        i__2 = l + l * h_dim1;
        h11.real = h__[i__2].real;
        h11.imag = h__[i__2].imag; // , expr subst
        i__2 = l + 1 + (l + 1) * h_dim1;
        h22.real = h__[i__2].real;
        h22.imag = h__[i__2].imag; // , expr subst
        z__1.real = h11.real - t.real;
        z__1.imag = h11.imag - t.imag; // , expr subst
        h11s.real = z__1.real;
        h11s.imag = z__1.imag; // , expr subst
        i__2 = l + 1 + l * h_dim1;
        h21 = h__[i__2].real;
        s = (d__1 = h11s.real, f2c_dabs(d__1)) + (d__2 = d_imag(&h11s), f2c_dabs(d__2))
            + f2c_dabs(h21);
        z__1.real = h11s.real / s;
        z__1.imag = h11s.imag / s; // , expr subst
        h11s.real = z__1.real;
        h11s.imag = z__1.imag; // , expr subst
        h21 /= s;
        v[0].real = h11s.real;
        v[0].imag = h11s.imag; // , expr subst
        v[1].real = h21;
        v[1].imag = 0.; // , expr subst
    L70: /* Single-shift QR step */
        i__2 = i__ - 1;
        for(k = m; k <= i__2; ++k)
        {
            /* The first iteration of this loop determines a reflection G */
            /* from the vector V and applies it from left and right to H, */
            /* thus creating a nonzero bulge below the subdiagonal. */
            /* Each subsequent iteration determines a reflection G to */
            /* restore the Hessenberg form in the (K-1)th column, and thus */
            /* chases the bulge one step toward the bottom of the active */
            /* submatrix. */
            /* V(2) is always real before the call to ZLARFG, and hence */
            /* after the call T2 ( = T1*V(2) ) is also real. */
            if(k > m)
            {
                aocl_blas_zcopy(&c__2, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
            }
            aocl_lapack_zlarfg(&c__2, v, &v[1], &c__1, &t1);
            if(k > m)
            {
                i__3 = k + (k - 1) * h_dim1;
                h__[i__3].real = v[0].real;
                h__[i__3].imag = v[0].imag; // , expr subst
                i__3 = k + 1 + (k - 1) * h_dim1;
                h__[i__3].real = 0.;
                h__[i__3].imag = 0.; // , expr subst
            }
            v2.real = v[1].real;
            v2.imag = v[1].imag; // , expr subst
            z__1.real = t1.real * v2.real - t1.imag * v2.imag;
            z__1.imag = t1.real * v2.imag + t1.imag * v2.real; // , expr subst
            t2 = z__1.real;
            /* Apply G from the left to transform the rows of the matrix */
            /* in columns K to I2. */
            i__3 = i2;
            for(j = k; j <= i__3; ++j)
            {
                d_cnjg(&z__3, &t1);
                i__4 = k + j * h_dim1;
                z__2.real = z__3.real * h__[i__4].real - z__3.imag * h__[i__4].imag;
                z__2.imag = z__3.real * h__[i__4].imag + z__3.imag * h__[i__4].real; // , expr subst
                i__5 = k + 1 + j * h_dim1;
                z__4.real = t2 * h__[i__5].real;
                z__4.imag = t2 * h__[i__5].imag; // , expr subst
                z__1.real = z__2.real + z__4.real;
                z__1.imag = z__2.imag + z__4.imag; // , expr subst
                sum.real = z__1.real;
                sum.imag = z__1.imag; // , expr subst
                i__4 = k + j * h_dim1;
                i__5 = k + j * h_dim1;
                z__1.real = h__[i__5].real - sum.real;
                z__1.imag = h__[i__5].imag - sum.imag; // , expr subst
                h__[i__4].real = z__1.real;
                h__[i__4].imag = z__1.imag; // , expr subst
                i__4 = k + 1 + j * h_dim1;
                i__5 = k + 1 + j * h_dim1;
                z__2.real = sum.real * v2.real - sum.imag * v2.imag;
                z__2.imag = sum.real * v2.imag + sum.imag * v2.real; // , expr subst
                z__1.real = h__[i__5].real - z__2.real;
                z__1.imag = h__[i__5].imag - z__2.imag; // , expr subst
                h__[i__4].real = z__1.real;
                h__[i__4].imag = z__1.imag; // , expr subst
                /* L80: */
            }
            /* Apply G from the right to transform the columns of the */
            /* matrix in rows I1 to fla_min(K+2,I). */
            /* Computing MIN */
            i__4 = k + 2;
            i__3 = fla_min(i__4, i__);
            for(j = i1; j <= i__3; ++j)
            {
                i__4 = j + k * h_dim1;
                z__2.real = t1.real * h__[i__4].real - t1.imag * h__[i__4].imag;
                z__2.imag = t1.real * h__[i__4].imag + t1.imag * h__[i__4].real; // , expr subst
                i__5 = j + (k + 1) * h_dim1;
                z__3.real = t2 * h__[i__5].real;
                z__3.imag = t2 * h__[i__5].imag; // , expr subst
                z__1.real = z__2.real + z__3.real;
                z__1.imag = z__2.imag + z__3.imag; // , expr subst
                sum.real = z__1.real;
                sum.imag = z__1.imag; // , expr subst
                i__4 = j + k * h_dim1;
                i__5 = j + k * h_dim1;
                z__1.real = h__[i__5].real - sum.real;
                z__1.imag = h__[i__5].imag - sum.imag; // , expr subst
                h__[i__4].real = z__1.real;
                h__[i__4].imag = z__1.imag; // , expr subst
                i__4 = j + (k + 1) * h_dim1;
                i__5 = j + (k + 1) * h_dim1;
                d_cnjg(&z__3, &v2);
                z__2.real = sum.real * z__3.real - sum.imag * z__3.imag;
                z__2.imag = sum.real * z__3.imag + sum.imag * z__3.real; // , expr subst
                z__1.real = h__[i__5].real - z__2.real;
                z__1.imag = h__[i__5].imag - z__2.imag; // , expr subst
                h__[i__4].real = z__1.real;
                h__[i__4].imag = z__1.imag; // , expr subst
                /* L90: */
            }
            if(*wantz)
            {
                /* Accumulate transformations in the matrix Z */
                i__3 = *ihiz;
                for(j = *iloz; j <= i__3; ++j)
                {
                    i__4 = j + k * z_dim1;
                    z__2.real = t1.real * z__[i__4].real - t1.imag * z__[i__4].imag;
                    z__2.imag = t1.real * z__[i__4].imag + t1.imag * z__[i__4].real; // , expr subst
                    i__5 = j + (k + 1) * z_dim1;
                    z__3.real = t2 * z__[i__5].real;
                    z__3.imag = t2 * z__[i__5].imag; // , expr subst
                    z__1.real = z__2.real + z__3.real;
                    z__1.imag = z__2.imag + z__3.imag; // , expr subst
                    sum.real = z__1.real;
                    sum.imag = z__1.imag; // , expr subst
                    i__4 = j + k * z_dim1;
                    i__5 = j + k * z_dim1;
                    z__1.real = z__[i__5].real - sum.real;
                    z__1.imag = z__[i__5].imag - sum.imag; // , expr subst
                    z__[i__4].real = z__1.real;
                    z__[i__4].imag = z__1.imag; // , expr subst
                    i__4 = j + (k + 1) * z_dim1;
                    i__5 = j + (k + 1) * z_dim1;
                    d_cnjg(&z__3, &v2);
                    z__2.real = sum.real * z__3.real - sum.imag * z__3.imag;
                    z__2.imag = sum.real * z__3.imag + sum.imag * z__3.real; // , expr subst
                    z__1.real = z__[i__5].real - z__2.real;
                    z__1.imag = z__[i__5].imag - z__2.imag; // , expr subst
                    z__[i__4].real = z__1.real;
                    z__[i__4].imag = z__1.imag; // , expr subst
                    /* L100: */
                }
            }
            if(k == m && m > l)
            {
                /* If the QR step was started at row M > L because two */
                /* consecutive small subdiagonals were found, then extra */
                /* scaling must be performed to ensure that H(M,M-1) remains */
                /* real. */
                z__1.real = 1. - t1.real;
                z__1.imag = 0. - t1.imag; // , expr subst
                temp.real = z__1.real;
                temp.imag = z__1.imag; // , expr subst
                d__1 = z_abs(&temp);
                z__1.real = temp.real / d__1;
                z__1.imag = temp.imag / d__1; // , expr subst
                temp.real = z__1.real;
                temp.imag = z__1.imag; // , expr subst
                i__3 = m + 1 + m * h_dim1;
                i__4 = m + 1 + m * h_dim1;
                d_cnjg(&z__2, &temp);
                z__1.real = h__[i__4].real * z__2.real - h__[i__4].imag * z__2.imag;
                z__1.imag = h__[i__4].real * z__2.imag + h__[i__4].imag * z__2.real; // , expr subst
                h__[i__3].real = z__1.real;
                h__[i__3].imag = z__1.imag; // , expr subst
                if(m + 2 <= i__)
                {
                    i__3 = m + 2 + (m + 1) * h_dim1;
                    i__4 = m + 2 + (m + 1) * h_dim1;
                    z__1.real = h__[i__4].real * temp.real - h__[i__4].imag * temp.imag;
                    z__1.imag = h__[i__4].real * temp.imag + h__[i__4].imag * temp.real; // , expr subst
                    h__[i__3].real = z__1.real;
                    h__[i__3].imag = z__1.imag; // , expr subst
                }
                i__3 = i__;
                for(j = m; j <= i__3; ++j)
                {
                    if(j != m + 1)
                    {
                        if(i2 > j)
                        {
                            i__4 = i2 - j;
                            aocl_blas_zscal(&i__4, &temp, &h__[j + (j + 1) * h_dim1], ldh);
                        }
                        i__4 = j - i1;
                        d_cnjg(&z__1, &temp);
                        aocl_blas_zscal(&i__4, &z__1, &h__[i1 + j * h_dim1], &c__1);
                        if(*wantz)
                        {
                            d_cnjg(&z__1, &temp);
                            aocl_blas_zscal(&nz, &z__1, &z__[*iloz + j * z_dim1], &c__1);
                        }
                    }
                    /* L110: */
                }
            }
            /* L120: */
        }
        /* Ensure that H(I,I-1) is real. */
        i__2 = i__ + (i__ - 1) * h_dim1;
        temp.real = h__[i__2].real;
        temp.imag = h__[i__2].imag; // , expr subst
        if(d_imag(&temp) != 0.)
        {
            rtemp = z_abs(&temp);
            i__2 = i__ + (i__ - 1) * h_dim1;
            h__[i__2].real = rtemp;
            h__[i__2].imag = 0.; // , expr subst
            z__1.real = temp.real / rtemp;
            z__1.imag = temp.imag / rtemp; // , expr subst
            temp.real = z__1.real;
            temp.imag = z__1.imag; // , expr subst
            if(i2 > i__)
            {
                i__2 = i2 - i__;
                d_cnjg(&z__1, &temp);
                aocl_blas_zscal(&i__2, &z__1, &h__[i__ + (i__ + 1) * h_dim1], ldh);
            }
            i__2 = i__ - i1;
            aocl_blas_zscal(&i__2, &temp, &h__[i1 + i__ * h_dim1], &c__1);
            if(*wantz)
            {
                aocl_blas_zscal(&nz, &temp, &z__[*iloz + i__ * z_dim1], &c__1);
            }
        }
        /* L130: */
    }
    /* Failure to converge in remaining number of iterations */
    *info = i__;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
L140: /* H(I,I-1) is negligible: one eigenvalue has converged. */
    i__1 = i__;
    i__2 = i__ + i__ * h_dim1;
    w[i__1].real = h__[i__2].real;
    w[i__1].imag = h__[i__2].imag; // , expr subst
    /* reset deflation counter */
    kdefl = 0;
    /* return to start of the main loop with new value of I. */
    i__ = l - 1;
    goto L30;
L150:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLAHQR */
}
/* zlahqr_ */
