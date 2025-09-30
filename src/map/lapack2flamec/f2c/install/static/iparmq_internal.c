#include "FLA_f2c.h" /* > \brief \b IPARMQ */


aocl_int64_t aocl_lapack_iparmq(aocl_int64_t *ispec, char *name__, char *opts, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                aocl_int64_t *lwork)
{
    /* System generated locals */
    aocl_int64_t ret_val, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double log(doublereal);
    /* Subroutine */
    int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Local variables */
    aocl_int64_t i__, ic, nh, ns, iz;
    char subnam[6];
    ftnlen name_len = strlen(name__);
    /* -- LAPACK auxiliary routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
    /* .. Scalar Arguments .. */
    /* ================================================================ */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    ns = 0;
    nh = 0;
    if(*ispec == 15 || *ispec == 13 || *ispec == 16)
    {
        /* ==== Set the number simultaneous shifts ==== */
        nh = *ihi - *ilo + 1;
        ns = 2;
        if(nh >= 30)
        {
            ns = 4;
        }
        if(nh >= 60)
        {
            ns = 10;
        }
        if(nh >= 150)
        {
            /* Computing MAX */
            r__1 = log((real)nh) / log(2.f);
            i__1 = 10;
            i__2 = nh / fla_i_nint(&r__1); // , expr subst
            ns = fla_max(i__1, i__2);
        }
        if(nh >= 590)
        {
            ns = 64;
        }
        if(nh >= 3000)
        {
            ns = 128;
        }
        if(nh >= 6000)
        {
            ns = 256;
        }
        /* Computing MAX */
        i__1 = 2;
        i__2 = ns - ns % 2; // , expr subst
        ns = fla_max(i__1, i__2);
    }
    if(*ispec == 12)
    {
        /* ===== Matrices of order smaller than NMIN get sent */
        /* . to xLAHQR, the classic double shift algorithm. */
        /* . This must be at least 11. ==== */
        ret_val = 75;
    }
    else if(*ispec == 14)
    {
        /* ==== INIBL: skip a multi-shift qr iteration and */
        /* . whenever aggressive early deflation finds */
        /* . at least (NIBBLE*(window size)/100) deflations. ==== */
        ret_val = 14;
    }
    else if(*ispec == 15)
    {
        /* ==== NSHFTS: The number of simultaneous shifts ===== */
        ret_val = ns;
    }
    else if(*ispec == 13)
    {
        /* ==== NW: deflation window size. ==== */
        if(nh <= 500)
        {
            ret_val = ns;
        }
        else
        {
            ret_val = ns * 3 / 2;
        }
    }
    else if(*ispec == 16)
    {
        /* ==== IACC22: Whether to accumulate reflections */
        /* . before updating the far-from-diagonal elements */
        /* . and whether to use 2-by-2 block structure while */
        /* . doing it. A small amount of work could be saved */
        /* . by making this choice dependent also upon the */
        /* . NH=IHI-ILO+1. */
        /* Convert NAME to upper case if the first character is lower case. */
        ret_val = 0;
        s_copy(subnam, name__, 6, name_len);
        ic = *(unsigned char *)subnam;
        iz = 'Z';
        if(iz == 90 || iz == 122)
        {
            /* ASCII character set */
            if(ic >= 97 && ic <= 122)
            {
                *(unsigned char *)subnam = (char)(ic - 32);
                for(i__ = 2; i__ <= 6; ++i__)
                {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if(ic >= 97 && ic <= 122)
                    {
                        *(unsigned char *)&subnam[i__ - 1] = (char)(ic - 32);
                    }
                }
            }
        }
        else if(iz == 233 || iz == 169)
        {
            /* EBCDIC character set */
            if(ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && ic <= 169)
            {
                *(unsigned char *)subnam = (char)(ic + 64);
                for(i__ = 2; i__ <= 6; ++i__)
                {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if(ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && ic <= 169)
                    {
                        *(unsigned char *)&subnam[i__ - 1] = (char)(ic + 64);
                    }
                }
            }
        }
        else if(iz == 218 || iz == 250)
        {
            /* Prime machines: ASCII+128 */
            if(ic >= 225 && ic <= 250)
            {
                *(unsigned char *)subnam = (char)(ic - 32);
                for(i__ = 2; i__ <= 6; ++i__)
                {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if(ic >= 225 && ic <= 250)
                    {
                        *(unsigned char *)&subnam[i__ - 1] = (char)(ic - 32);
                    }
                }
            }
        }
        if(s_cmp(subnam + 1, "GGHRD", (ftnlen)5, (ftnlen)5) == 0
           || s_cmp(subnam + 1, "GGHD3", (ftnlen)5, (ftnlen)5) == 0)
        {
            ret_val = 1;
            if(nh >= 14)
            {
                ret_val = 2;
            }
        }
        else if(s_cmp(subnam + 3, "EXC", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(nh >= 14)
            {
                ret_val = 1;
            }
            if(nh >= 14)
            {
                ret_val = 2;
            }
        }
        else if(s_cmp(subnam + 1, "HSEQR", (ftnlen)5, (ftnlen)5) == 0
                || s_cmp(subnam + 1, "LAQR", (ftnlen)4, (ftnlen)4) == 0)
        {
            if(ns >= 14)
            {
                ret_val = 1;
            }
            if(ns >= 14)
            {
                ret_val = 2;
            }
        }
    }
    else
    {
        /* ===== invalid value of ispec ===== */
        ret_val = -1;
    }
    /* ==== End of IPARMQ ==== */
    return ret_val;
}
/* iparmq_ */
