#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static integer c__1 = 1;
static integer c_n1 = -1;

int dgesdd_check(char *jobz, integer *m, integer *n, double *a, integer *lda, double *s, double *u,
                 integer *ldu, double *vt, integer *ldvt, double *work, integer *lwork,
                 integer *iwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    /* Local variables */
    integer minmn, wrkbl, mnthr;
    logical wntqa;
    logical wntqn, wntqo, wntqs;
    integer bdspac;
    integer minwrk, maxwrk;
    logical wntqas, lquery;
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "dgesdd inputs: jobz %c, m %d, n %d, lda %d, ldu %d, ldvt %d\n", *jobz, *m, *n,
            *lda, *ldu, *ldvt);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    minmn = fla_min(*m, *n);
    wntqa = lsame_(jobz, "A", 1, 1);
    wntqs = lsame_(jobz, "S", 1, 1);
    wntqas = wntqa || wntqs;
    wntqo = lsame_(jobz, "O", 1, 1);
    wntqn = lsame_(jobz, "N", 1, 1);
    lquery = *lwork == -1;
    if(!(wntqa || wntqs || wntqo || wntqn))
    {
        *info = -1;
    }
    else if(*m < 0)
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    else if(*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *m)
    {
        *info = -8;
    }
    else if(*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn
            || wntqo && *m >= *n && *ldvt < *n)
    {
        *info = -10;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV.) */
    if(*info == 0)
    {
        minwrk = 1;
        maxwrk = 1;
        if(*m >= *n && minmn > 0)
        {
            /* Compute space needed for DBDSDC */
            mnthr = (integer)(minmn * 11. / 6.);
            if(wntqn)
            {
                bdspac = *n * 7;
            }
            else
            {
                bdspac = *n * 3 * *n + (*n << 2);
            }
            if(*m >= mnthr)
            {
                if(wntqn)
                {
                    /* Path 1 (M much larger than N, JOBZ='N') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *n * 3
                          + (*n << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = bdspac + *n;
                }
                else if(wntqo)
                {
                    /* Path 2 (M much larger than N, JOBZ='O') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *n * ilaenv_(&c__1, "DORGQR", " ", m, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *n * 3
                          + (*n << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + (*n << 1) * *n;
                    minwrk = bdspac + (*n << 1) * *n + *n * 3;
                }
                else if(wntqs)
                {
                    /* Path 3 (M much larger than N, JOBZ='S') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *n * ilaenv_(&c__1, "DORGQR", " ", m, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *n * 3
                          + (*n << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + *n * *n;
                    minwrk = bdspac + *n * *n + *n * 3;
                }
                else if(wntqa)
                {
                    /* Path 4 (M much larger than N, JOBZ='A') */
                    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n + *m * ilaenv_(&c__1, "DORGQR", " ", m, m, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *n * 3
                          + (*n << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "QLN", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + *n * *n;
                    minwrk = bdspac + *n * *n + (*n << 1) + *m;
                }
            }
            else
            {
                /* Path 5 (M at least N, but not much larger) */
                wrkbl = *n * 3 + (*m + *n) * ilaenv_(&c__1, "DGEBRD", " ", m, n, &c_n1, &c_n1);
                if(wntqn)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = *n * 3 + fla_max(*m, bdspac);
                }
                else if(wntqo)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "QLN", m, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + *m * *n;
                    /* Computing MAX */
                    i__1 = *m;
                    i__2 = *n * *n + bdspac; // , expr subst
                    minwrk = *n * 3 + fla_max(i__1, i__2);
                }
                else if(wntqs)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "QLN", m, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *n * 3; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = *n * 3 + fla_max(*m, bdspac);
                }
                else if(wntqa)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *n * 3
                           + *n * ilaenv_(&c__1, "DORMBR", "PRT", n, n, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = bdspac + *n * 3; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = *n * 3 + fla_max(*m, bdspac);
                }
            }
        }
        else if(minmn > 0)
        {
            /* Compute space needed for DBDSDC */
            mnthr = (integer)(minmn * 11. / 6.);
            if(wntqn)
            {
                bdspac = *m * 7;
            }
            else
            {
                bdspac = *m * 3 * *m + (*m << 2);
            }
            if(*n >= mnthr)
            {
                if(wntqn)
                {
                    /* Path 1t (N much larger than M, JOBZ='N') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *m * 3
                          + (*m << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = bdspac + *m;
                }
                else if(wntqo)
                {
                    /* Path 2t (N much larger than M, JOBZ='O') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *m * ilaenv_(&c__1, "DORGLQ", " ", m, n, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *m * 3
                          + (*m << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "PRT", m, m, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + (*m << 1) * *m;
                    minwrk = bdspac + (*m << 1) * *m + *m * 3;
                }
                else if(wntqs)
                {
                    /* Path 3t (N much larger than M, JOBZ='S') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *m * ilaenv_(&c__1, "DORGLQ", " ", m, n, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *m * 3
                          + (*m << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "PRT", m, m, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + *m * *m;
                    minwrk = bdspac + *m * *m + *m * 3;
                }
                else if(wntqa)
                {
                    /* Path 4t (N much larger than M, JOBZ='A') */
                    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m + *n * ilaenv_(&c__1, "DORGLQ", " ", n, n, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2
                        = *m * 3
                          + (*m << 1)
                                * ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "QLN", m, m, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "PRT", m, m, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + *m * *m;
                    minwrk = bdspac + *m * *m + *m * 3;
                }
            }
            else
            {
                /* Path 5t (N greater than M, but not much larger) */
                wrkbl = *m * 3 + (*m + *n) * ilaenv_(&c__1, "DGEBRD", " ", m, n, &c_n1, &c_n1);
                if(wntqn)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = *m * 3 + fla_max(*n, bdspac);
                }
                else if(wntqo)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "PRT", m, n, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    maxwrk = wrkbl + *m * *n;
                    /* Computing MAX */
                    i__1 = *n;
                    i__2 = *m * *m + bdspac; // , expr subst
                    minwrk = *m * 3 + fla_max(i__1, i__2);
                }
                else if(wntqs)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "PRT", m, n, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = *m * 3 + fla_max(*n, bdspac);
                }
                else if(wntqa)
                {
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "QLN", m, m, n, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = *m * 3
                           + *m * ilaenv_(&c__1, "DORMBR", "PRT", n, n, m, &c_n1); // , expr subst
                    wrkbl = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = wrkbl;
                    i__2 = bdspac + *m * 3; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    minwrk = *m * 3 + fla_max(*n, bdspac);
                }
            }
        }
        maxwrk = fla_max(maxwrk, minwrk);
        work[1] = (double)maxwrk;
        if(*lwork < minwrk && !lquery)
        {
            *info = -12;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGESDD", &i__1, (ftnlen)6);
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
