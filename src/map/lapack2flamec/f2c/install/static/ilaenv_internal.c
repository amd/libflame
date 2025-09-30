#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static real c_b179 = 0.f;
static real c_b180 = 1.f;
static aocl_int64_t c__0 = 0;

aocl_int64_t aocl_lapack_ilaenv(aocl_int64_t *ispec, char *name__, char *opts, aocl_int64_t *n1, aocl_int64_t *n2, aocl_int64_t *n3,
                aocl_int64_t *n4)
{
    /* System generated locals */
    aocl_int64_t ret_val, i__1, i__2, i__3;
    /* Builtin functions */
    /* Subroutine */
    int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_len(char *), s_cmp(char *, char *, ftnlen, ftnlen);
    /* strlen(s1)* Local variables */
    logical twostage;
    aocl_int64_t i__;
    char c1[1], c2[2], c3[3], c4[2];
    aocl_int64_t ic, nb, iz, nx;
    logical cname;
    aocl_int64_t nbmin;
    logical sname;
    char subnam[16];
    ftnlen name_len = strlen(name__);
    /* -- LAPACK auxiliary routine (version 3.9.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
     */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    switch(*ispec)
    {
        case 1:
            goto L10;
        case 2:
            goto L10;
        case 3:
            goto L10;
        case 4:
            goto L80;
        case 5:
            goto L90;
        case 6:
            goto L100;
        case 7:
            goto L110;
        case 8:
            goto L120;
        case 9:
            goto L130;
        case 10:
            goto L140;
        case 11:
            goto L150;
        case 12:
            goto L160;
        case 13:
            goto L160;
        case 14:
            goto L160;
        case 15:
            goto L160;
        case 16:
            goto L160;
        case 17:
            goto L160;
    }
    /* Invalid value for ISPEC */
    ret_val = -1;
    return ret_val;
L10: /* Convert NAME to upper case if the first character is lower case. */
    ret_val = 1;
    s_copy(subnam, name__, (ftnlen)16, name_len);
    /*if(name_len < 15)
      subnam[name_len] = '\0';*/
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
                /* L20: */
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
                /* L30: */
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
                /* L40: */
            }
        }
    }
    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if(!(cname || sname))
    {
        return ret_val;
    }
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);
    twostage = i_len(subnam) >= 11 && *(unsigned char *)&subnam[10] == '2';
    switch(*ispec)
    {
        case 1:
            goto L50;
        case 2:
            goto L60;
        case 3:
            goto L70;
    }
L50: /* ISPEC = 1: block size */
    /* In these examples, separate code is provided for setting NB for */
    /* real and scomplex. We assume that NB will take the same value in */
    /* single or double precision. */
    nb = 1;
    if(s_cmp(subnam + 1, "LAORH", (ftnlen)5, (ftnlen)5) == 0)
    {
        /* This is for *LAORHR_GETRFNP routine */
        if(sname)
        {
            nb = 32;
        }
        else
        {
            nb = 32;
        }
    }
    else if(s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 64;
            }
            else
            {
                nb = 64;
            }
        }
        else if(s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0
                || s_cmp(c3, "RQF", (ftnlen)3, (ftnlen)3) == 0
                || s_cmp(c3, "LQF", (ftnlen)3, (ftnlen)3) == 0
                || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 32;
            }
            else
            {
                nb = 32;
            }
        }
        else if(s_cmp(c3, "QR", (ftnlen)2, (ftnlen)2) == 0)
        {
            if(*n3 == 1)
            {
                if(sname)
                {
                    /* M*N */
                    if(*n1 * *n2 <= 131072 || *n1 <= 8192)
                    {
                        nb = *n1;
                    }
                    else
                    {
                        nb = 32768 / *n2;
                    }
                }
                else
                {
                    if(*n1 * *n2 <= 131072 || *n1 <= 8192)
                    {
                        nb = *n1;
                    }
                    else
                    {
                        nb = 32768 / *n2;
                    }
                }
            }
            else
            {
                if(sname)
                {
                    nb = 1;
                }
                else
                {
                    nb = 1;
                }
            }
        }
        else if(s_cmp(c3, "LQ", (ftnlen)2, (ftnlen)2) == 0)
        {
            if(*n3 == 2)
            {
                if(sname)
                {
                    /* M*N */
                    if(*n1 * *n2 <= 131072 || *n1 <= 8192)
                    {
                        nb = *n1;
                    }
                    else
                    {
                        nb = 32768 / *n2;
                    }
                }
                else
                {
                    if(*n1 * *n2 <= 131072 || *n1 <= 8192)
                    {
                        nb = *n1;
                    }
                    else
                    {
                        nb = 32768 / *n2;
                    }
                }
            }
            else
            {
                if(sname)
                {
                    nb = 1;
                }
                else
                {
                    nb = 1;
                }
            }
        }
        else if(s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 32;
            }
            else
            {
                nb = 32;
            }
        }
        else if(s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 32;
            }
            else
            {
                nb = 32;
            }
        }
        else if(s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 64;
            }
            else
            {
                nb = 64;
            }
        }
        else if(s_cmp(subnam + 3, "QP3RK", (ftnlen)5, (ftnlen)5) == 0)
        {
            if(sname)
            {
                nb = 32;
            }
            else
            {
                nb = 32;
            }
        }
    }
    else if(s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 64;
            }
            else
            {
                nb = 64;
            }
        }
    }
    else if(s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                if(twostage)
                {
                    nb = 192;
                }
                else
                {
                    nb = 64;
                }
            }
            else
            {
                if(twostage)
                {
                    nb = 192;
                }
                else
                {
                    nb = 64;
                }
            }
        }
        else if(sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            nb = 32;
        }
        else if(sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0)
        {
            nb = 64;
        }
    }
    else if(cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(twostage)
            {
                nb = 192;
            }
            else
            {
                nb = 64;
            }
        }
        else if(s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            nb = 32;
        }
        else if(s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0)
        {
            nb = 64;
        }
    }
    else if(sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(*(unsigned char *)c3 == 'G')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nb = 32;
            }
        }
        else if(*(unsigned char *)c3 == 'M')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nb = 32;
            }
        }
    }
    else if(cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(*(unsigned char *)c3 == 'G')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nb = 32;
            }
        }
        else if(*(unsigned char *)c3 == 'M')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nb = 32;
            }
        }
    }
    else if(s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                if(*n4 <= 64)
                {
                    nb = 1;
                }
                else
                {
                    nb = 32;
                }
            }
            else
            {
                if(*n4 <= 64)
                {
                    nb = 1;
                }
                else
                {
                    nb = 32;
                }
            }
        }
    }
    else if(s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                if(*n2 <= 64)
                {
                    nb = 1;
                }
                else
                {
                    nb = 32;
                }
            }
            else
            {
                if(*n2 <= 64)
                {
                    nb = 1;
                }
                else
                {
                    nb = 32;
                }
            }
        }
    }
    else if(s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 64;
            }
            else
            {
                nb = 64;
            }
        }
        else if(s_cmp(c3, "EVC", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 64;
            }
            else
            {
                nb = 64;
            }
        }
        else if(s_cmp(c3, "SYL", (ftnlen)3, (ftnlen)3) == 0)
        {
            /* The upper bound is to prevent overly aggressive scaling. */
            if(sname)
            {
                /* Computing MIN */
                /* Computing MAX */
                i__2 = 48;
                i__3 = (fla_min(*n1, *n2) << 4) / 100; // , expr subst
                i__1 = fla_max(i__2, i__3);
                nb = fla_min(i__1, 240);
            }
            else
            {
                /* Computing MIN */
                /* Computing MAX */
                i__2 = 24;
                i__3 = (fla_min(*n1, *n2) << 3) / 100; // , expr subst
                i__1 = fla_max(i__2, i__3);
                nb = fla_min(i__1, 80);
            }
        }
    }
    else if(s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 64;
            }
            else
            {
                nb = 64;
            }
        }
        else if(s_cmp(c3, "TRS", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 32;
            }
            else
            {
                nb = 32;
            }
        }
    }
    else if(sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0)
        {
            nb = 1;
        }
    }
    else if(s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0)
    {
        nb = 32;
        if(s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nb = 32;
            }
            else
            {
                nb = 32;
            }
        }
    }
    ret_val = nb;
    return ret_val;
L60: /* ISPEC = 2: minimum block size */
    nbmin = 2;
    if(s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0
           || s_cmp(c3, "RQF", (ftnlen)3, (ftnlen)3) == 0
           || s_cmp(c3, "LQF", (ftnlen)3, (ftnlen)3) == 0
           || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nbmin = 2;
            }
            else
            {
                nbmin = 2;
            }
        }
        else if(s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nbmin = 2;
            }
            else
            {
                nbmin = 2;
            }
        }
        else if(s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nbmin = 2;
            }
            else
            {
                nbmin = 2;
            }
        }
        else if(s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nbmin = 2;
            }
            else
            {
                nbmin = 2;
            }
        }
        else if(s_cmp(subnam + 3, "QP3RK", (ftnlen)5, (ftnlen)5) == 0)
        {
            if(sname)
            {
                nbmin = 2;
            }
            else
            {
                nbmin = 2;
            }
        }
    }
    else if(s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nbmin = 8;
            }
            else
            {
                nbmin = 8;
            }
        }
        else if(sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            nbmin = 2;
        }
    }
    else if(cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            nbmin = 2;
        }
    }
    else if(sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(*(unsigned char *)c3 == 'G')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nbmin = 2;
            }
        }
        else if(*(unsigned char *)c3 == 'M')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nbmin = 2;
            }
        }
    }
    else if(cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(*(unsigned char *)c3 == 'G')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nbmin = 2;
            }
        }
        else if(*(unsigned char *)c3 == 'M')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nbmin = 2;
            }
        }
    }
    else if(s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0)
    {
        nbmin = 2;
        if(s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0)
        {
            nbmin = 2;
        }
    }
    ret_val = nbmin;
    return ret_val;
L70: /* ISPEC = 3: crossover point */
    nx = 0;
    if(s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0
           || s_cmp(c3, "RQF", (ftnlen)3, (ftnlen)3) == 0
           || s_cmp(c3, "LQF", (ftnlen)3, (ftnlen)3) == 0
           || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nx = 128;
            }
            else
            {
                nx = 128;
            }
        }
        else if(s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nx = 128;
            }
            else
            {
                nx = 128;
            }
        }
        else if(s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            if(sname)
            {
                nx = 128;
            }
            else
            {
                nx = 128;
            }
        }
        else if(s_cmp(subnam + 3, "QP3RK", (ftnlen)5, (ftnlen)5) == 0)
        {
            if(sname)
            {
                nx = 128;
            }
            else
            {
                nx = 128;
            }
        }
    }
    else if(s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            nx = 32;
        }
    }
    else if(cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0)
        {
            nx = 32;
        }
    }
    else if(sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(*(unsigned char *)c3 == 'G')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nx = 128;
            }
        }
    }
    else if(cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0)
    {
        if(*(unsigned char *)c3 == 'G')
        {
            if(s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "RQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "LQ", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "TR", (ftnlen)2, (ftnlen)2) == 0
               || s_cmp(c4, "BR", (ftnlen)2, (ftnlen)2) == 0)
            {
                nx = 128;
            }
        }
    }
    else if(s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0)
    {
        nx = 128;
        if(s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0)
        {
            nx = 128;
        }
    }
    ret_val = nx;
    return ret_val;
L80: /* ISPEC = 4: number of shifts (used by xHSEQR) */
    ret_val = 6;
    return ret_val;
L90: /* ISPEC = 5: minimum column dimension (not used) */
    ret_val = 2;
    return ret_val;
L100: /* ISPEC = 6: crossover point for SVD (used by xGELSS and xGESVD) */
    ret_val = (integer)((real)fla_min(*n1, *n2) * 1.6f);
    return ret_val;
L110: /* ISPEC = 7: number of processors (not used) */
    ret_val = 1;
    return ret_val;
L120: /* ISPEC = 8: crossover point for multishift (used by xHSEQR) */
    ret_val = 50;
    return ret_val;
L130: /* ISPEC = 9: maximum size of the subproblems at the bottom of the */
    /* computation tree in the divide-and-conquer algorithm */
    /* (used by xGELSD and xGESDD) */
    ret_val = 25;
    return ret_val;
L140: /* ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
       */
    /* ILAENV = 0 */
    ret_val = 1;
    if(ret_val == 1)
    {
        ret_val = aocl_lapack_ieeeck(&c__1, &c_b179, &c_b180);
    }
    return ret_val;
L150: /* ISPEC = 11: ieee infinity arithmetic can be trusted not to trap */
    /* ILAENV = 0 */
    ret_val = 1;
    if(ret_val == 1)
    {
        ret_val = aocl_lapack_ieeeck(&c__0, &c_b179, &c_b180);
    }
    return ret_val;
L160: /* 12 <= ISPEC <= 16: xHSEQR or related subroutines. */
    ret_val = aocl_lapack_iparmq(ispec, name__, opts, n1, n2, n3, n4); // , name_len, opts_len) ;
    return ret_val;
    /* End of ILAENV */
}
/* ilaenv_ */
