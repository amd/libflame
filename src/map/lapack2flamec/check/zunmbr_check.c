#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static integer c__1 = 1;
static integer c_n1 = -1;

int zunmbr_check(char *vect, char *side, char *trans, integer *m, integer *n, integer *k,
                 dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc,
                 dcomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    char ch__1[2];

    /* Local variables */
    integer nb, nq, nw;
    logical left;
    logical notran, applyq;
    integer lwkopt;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    applyq = lsame_(vect, "Q", 1, 1);
    left = lsame_(side, "L", 1, 1);
    notran = lsame_(trans, "N", 1, 1);
    lquery = *lwork == -1;
    /* NQ is the order of Q or P and NW is the minimum dimension of WORK */
    if(left)
    {
        nq = *m;
        nw = *n;
    }
    else
    {
        nq = *n;
        nw = *m;
    }
    if(*m == 0 || *n == 0)
    {
        nw = 0;
    }
    if(!applyq && !lsame_(vect, "P", 1, 1))
    {
        *info = -1;
    }
    else if(!left && !lsame_(side, "R", 1, 1))
    {
        *info = -2;
    }
    else if(!notran && !lsame_(trans, "C", 1, 1))
    {
        *info = -3;
    }
    else if(*m < 0)
    {
        *info = -4;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*k < 0)
    {
        *info = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = fla_min(nq, *k); // , expr subst
        if(applyq && *lda < fla_max(1, nq) || !applyq && *lda < fla_max(i__1, i__2))
        {
            *info = -8;
        }
        else if(*ldc < fla_max(1, *m))
        {
            *info = -11;
        }
        else if(*lwork < fla_max(1, nw) && !lquery)
        {
            *info = -13;
        }
    }
    if(*info == 0)
    {
        if(nw > 0)
        {
            if(applyq)
            {
                if(left)
                {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    nb = ilaenv_(&c__1, "ZUNMQR", ch__1, &i__1, n, &i__2, &c_n1);
                }
                else
                {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    nb = ilaenv_(&c__1, "ZUNMQR", ch__1, m, &i__1, &i__2, &c_n1);
                }
            }
            else
            {
                if(left)
                {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    nb = ilaenv_(&c__1, "ZUNMLQ", ch__1, &i__1, n, &i__2, &c_n1);
                }
                else
                {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    nb = ilaenv_(&c__1, "ZUNMLQ", ch__1, m, &i__1, &i__2, &c_n1);
                }
            }
            /* Computing MAX */
            i__1 = 1;
            i__2 = nw * nb; // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        else
        {
            lwkopt = 1;
        }
        work[1].real = (double)lwkopt;
        work[1].imag = 0.; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNMBR", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
