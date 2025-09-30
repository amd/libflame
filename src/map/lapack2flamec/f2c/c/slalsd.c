/* ../netlib/slalsd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static real c_b6 = 0.f;
static aocl_int64_t c__0 = 0;
static real c_b11 = 1.f;
/* > \brief \b SLALSD uses the singular value decomposition of A to solve the least squares problem.
 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLALSD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slalsd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slalsd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slalsd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, */
/* RANK, WORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDB, N, NRHS, RANK, SMLSIZ */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL B( LDB, * ), D( * ), E( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLALSD uses the singular value decomposition of A to solve the least */
/* > squares problem of finding X to minimize the Euclidean norm of each */
/* > column of A*X-B, where A is N-by-N upper bidiagonal, and X and B */
/* > are N-by-NRHS. The solution X overwrites B. */
/* > */
/* > The singular values of A smaller than RCOND times the largest */
/* > singular value are treated as zero in solving the least squares */
/* > problem;
in this case a minimum norm solution is returned. */
/* > The actual singular values are returned in D in ascending order. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': D and E define an upper bidiagonal matrix. */
/* > = 'L': D and E define a lower bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLSIZ */
/* > \verbatim */
/* > SMLSIZ is INTEGER */
/* > The maximum size of the subproblems at the bottom of the */
/* > computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the bidiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of columns of B. NRHS must be at least 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry D contains the main diagonal of the bidiagonal */
/* > matrix. On exit, if INFO = 0, D contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > Contains the super-diagonal entries of the bidiagonal matrix. */
/* > On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,NRHS) */
/* > On input, B contains the right hand sides of the least */
/* > squares problem. On output, B contains the solution X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B in the calling subprogram. */
/* > LDB must be at least fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The singular values of A less than or equal to RCOND times */
/* > the largest singular value are treated as zero in solving */
/* > the least squares problem. If RCOND is negative, */
/* > machine precision is used instead. */
/* > For example, if diag(S)*X=B were the least squares problem, */
/* > where diag(S) is a diagonal matrix of singular values, the */
/* > solution would be X(i) = B(i) / S(i) if S(i) is greater than */
/* > RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to */
/* > RCOND*max(S). */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* > RANK is INTEGER */
/* > The number of singular values of A greater than RCOND times */
/* > the largest singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension at least */
/* > (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2), */
/* > where NLVL = fla_max(0, INT(log_2 (N/(SMLSIZ+1))) + 1). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension at least */
/* > (3*N*NLVL + 11*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: The algorithm failed to compute a singular value while */
/* > working on the submatrix lying in rows and columns */
/* > INFO/(N+1) through MOD(INFO,N+1). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* > California at Berkeley, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void slalsd_(char *uplo, aocl_int_t *smlsiz, aocl_int_t *n, aocl_int_t *nrhs, real *d__, real *e,
             real *b, aocl_int_t *ldb, real *rcond, aocl_int_t *rank, real *work, aocl_int_t *iwork,
             aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_slalsd(uplo, smlsiz, n, nrhs, d__, e, b, ldb, rcond, rank, work, iwork, info);
#else
    aocl_int64_t smlsiz_64 = *smlsiz;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t rank_64 = *rank;
    aocl_int64_t info_64 = *info;

    aocl_lapack_slalsd(uplo, &smlsiz_64, &n_64, &nrhs_64, d__, e, b, &ldb_64, rcond, &rank_64, work,
                       iwork, &info_64);

    *rank = (aocl_int_t)rank_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_slalsd(char *uplo, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs,
                        real *d__, real *e, real *b, aocl_int64_t *ldb, real *rcond,
                        aocl_int64_t *rank, real *work, aocl_int_t *iwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slalsd inputs: uplo %c ,smlsiz %" FLA_IS ",n %" FLA_IS ",nrhs %" FLA_IS
                      ",ldb %" FLA_IS "",
                      *uplo, *smlsiz, *n, *nrhs, *ldb);
    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double log(doublereal), r_sign(real *, real *);
    /* Local variables */
    aocl_int64_t c__, i__, j, k;
    real r__;
    aocl_int64_t s, u, z__;
    real cs;
    aocl_int64_t bx;
    real sn;
    aocl_int64_t st, vt, nm1, st1;
    real eps;
    aocl_int64_t iwk;
    real tol;
    aocl_int64_t difl, difr;
    real rcnd;
    aocl_int64_t perm, nsub, nlvl, sqre, bxst;
    aocl_int64_t poles, sizei, nsize;
    aocl_int64_t nwork, icmpq1, icmpq2;
    extern real slamch_(char *);
    aocl_int64_t givcol;
    real orgnrm;
    aocl_int64_t givnum;
    aocl_int64_t givptr, smlszp;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -3;
    }
    else if(*nrhs < 1)
    {
        *info = -4;
    }
    else if(*ldb < 1 || *ldb < *n)
    {
        *info = -8;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SLALSD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    eps = slamch_("Epsilon");
    /* Set up the tolerance. */
    if(*rcond <= 0.f || *rcond >= 1.f)
    {
        rcnd = eps;
    }
    else
    {
        rcnd = *rcond;
    }
    *rank = 0;
    /* Quick return if possible. */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*n == 1)
    {
        if(d__[1] == 0.f)
        {
            aocl_lapack_slaset("A", &c__1, nrhs, &c_b6, &c_b6, &b[b_offset], ldb);
        }
        else
        {
            *rank = 1;
            aocl_lapack_slascl("G", &c__0, &c__0, &d__[1], &c_b11, &c__1, nrhs, &b[b_offset], ldb,
                               info);
            d__[1] = f2c_abs(d__[1]);
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Rotate the matrix if it is lower bidiagonal. */
    if(*(unsigned char *)uplo == 'L')
    {
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            slartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            if(*nrhs == 1)
            {
                aocl_blas_srot(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &c__1, &cs,
                               &sn);
            }
            else
            {
                work[(i__ << 1) - 1] = cs;
                work[i__ * 2] = sn;
            }
            /* L10: */
        }
        if(*nrhs > 1)
        {
            i__1 = *nrhs;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = *n - 1;
                for(j = 1; j <= i__2; ++j)
                {
                    cs = work[(j << 1) - 1];
                    sn = work[j * 2];
                    aocl_blas_srot(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ * b_dim1],
                                   &c__1, &cs, &sn);
                    /* L20: */
                }
                /* L30: */
            }
        }
    }
    /* Scale. */
    nm1 = *n - 1;
    orgnrm = aocl_lapack_slanst("M", n, &d__[1], &e[1]);
    if(orgnrm == 0.f)
    {
        aocl_lapack_slaset("A", n, nrhs, &c_b6, &c_b6, &b[b_offset], ldb);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    aocl_lapack_slascl("G", &c__0, &c__0, &orgnrm, &c_b11, n, &c__1, &d__[1], n, info);
    aocl_lapack_slascl("G", &c__0, &c__0, &orgnrm, &c_b11, &nm1, &c__1, &e[1], &nm1, info);
    /* If N is smaller than the minimum divide size SMLSIZ, then solve */
    /* the problem with another solver. */
    if(*n <= *smlsiz)
    {
        nwork = *n * *n + 1;
        aocl_lapack_slaset("A", n, n, &c_b6, &c_b11, &work[1], n);
        aocl_lapack_slasdq("U", &c__0, n, n, &c__0, nrhs, &d__[1], &e[1], &work[1], n, &work[1], n,
                           &b[b_offset], ldb, &work[nwork], info);
        if(*info != 0)
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        tol = rcnd * (r__1 = d__[aocl_blas_isamax(n, &d__[1], &c__1)], f2c_abs(r__1));
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if(d__[i__] <= tol)
            {
                aocl_lapack_slaset("A", &c__1, nrhs, &c_b6, &c_b6, &b[i__ + b_dim1], ldb);
            }
            else
            {
                aocl_lapack_slascl("G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs,
                                   &b[i__ + b_dim1], ldb, info);
                ++(*rank);
            }
            /* L40: */
        }
        aocl_blas_sgemm("T", "N", n, nrhs, n, &c_b11, &work[1], n, &b[b_offset], ldb, &c_b6,
                        &work[nwork], n);
        aocl_lapack_slacpy("A", n, nrhs, &work[nwork], n, &b[b_offset], ldb);
        /* Unscale. */
        aocl_lapack_slascl("G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, info);
        aocl_lapack_slasrt("D", n, &d__[1], info);
        aocl_lapack_slascl("G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], ldb, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Book-keeping and setting up some constants. */
    nlvl = (integer)(log((real)(*n) / (real)(*smlsiz + 1)) / log(2.f)) + 1;
    smlszp = *smlsiz + 1;
    u = 1;
    vt = *smlsiz * *n + 1;
    difl = vt + smlszp * *n;
    difr = difl + nlvl * *n;
    z__ = difr + (nlvl * *n << 1);
    c__ = z__ + nlvl * *n;
    s = c__ + *n;
    poles = s + *n;
    givnum = poles + (nlvl << 1) * *n;
    bx = givnum + (nlvl << 1) * *n;
    nwork = bx + *n * *nrhs;
    sizei = *n + 1;
    k = sizei + *n;
    givptr = k + *n;
    perm = givptr + *n;
    givcol = perm + nlvl * *n;
    iwk = givcol + (nlvl * *n << 1);
    st = 1;
    sqre = 0;
    icmpq1 = 1;
    icmpq2 = 0;
    nsub = 0;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if((r__1 = d__[i__], f2c_abs(r__1)) < eps)
        {
            d__[i__] = r_sign(&eps, &d__[i__]);
        }
        /* L50: */
    }
    i__1 = nm1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if((r__1 = e[i__], f2c_abs(r__1)) < eps || i__ == nm1)
        {
            ++nsub;
            iwork[nsub] = (aocl_int_t)(st);
            /* Subproblem found. First determine its size and then */
            /* apply divide and conquer on it. */
            if(i__ < nm1)
            {
                /* A subproblem with E(I) small for I < NM1. */
                nsize = i__ - st + 1;
                iwork[sizei + nsub - 1] = (aocl_int_t)(nsize);
            }
            else if((r__1 = e[i__], f2c_abs(r__1)) >= eps)
            {
                /* A subproblem with E(NM1) not too small but I = NM1. */
                nsize = *n - st + 1;
                iwork[sizei + nsub - 1] = (aocl_int_t)(nsize);
            }
            else
            {
                /* A subproblem with E(NM1) small. This implies an */
                /* 1-by-1 subproblem at D(N), which is not solved */
                /* explicitly. */
                nsize = i__ - st + 1;
                iwork[sizei + nsub - 1] = (aocl_int_t)(nsize);
                ++nsub;
                iwork[nsub] = (aocl_int_t)(*n);
                iwork[sizei + nsub - 1] = 1;
                aocl_blas_scopy(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
            }
            st1 = st - 1;
            if(nsize == 1)
            {
                /* This is a 1-by-1 subproblem and is not solved */
                /* explicitly. */
                aocl_blas_scopy(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
            }
            else if(nsize <= *smlsiz)
            {
                /* This is a small subproblem and is solved by SLASDQ. */
                aocl_lapack_slaset("A", &nsize, &nsize, &c_b6, &c_b11, &work[vt + st1], n);
                aocl_lapack_slasdq("U", &c__0, &nsize, &nsize, &c__0, nrhs, &d__[st], &e[st],
                                   &work[vt + st1], n, &work[nwork], n, &b[st + b_dim1], ldb,
                                   &work[nwork], info);
                if(*info != 0)
                {
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
                aocl_lapack_slacpy("A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
            }
            else
            {
                /* A large problem. Solve it using divide and conquer. */
                aocl_lapack_slasda(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &work[u + st1],
                                   n, &work[vt + st1], &iwork[k + st1], &work[difl + st1],
                                   &work[difr + st1], &work[z__ + st1], &work[poles + st1],
                                   &iwork[givptr + st1], &iwork[givcol + st1], n,
                                   &iwork[perm + st1], &work[givnum + st1], &work[c__ + st1],
                                   &work[s + st1], &work[nwork], &iwork[iwk], info);
                if(*info != 0)
                {
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
                bxst = bx + st1;
                aocl_lapack_slalsa(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &work[bxst],
                                   n, &work[u + st1], n, &work[vt + st1], &iwork[k + st1],
                                   &work[difl + st1], &work[difr + st1], &work[z__ + st1],
                                   &work[poles + st1], &iwork[givptr + st1], &iwork[givcol + st1],
                                   n, &iwork[perm + st1], &work[givnum + st1], &work[c__ + st1],
                                   &work[s + st1], &work[nwork], &iwork[iwk], info);
                if(*info != 0)
                {
                    AOCL_DTL_TRACE_LOG_EXIT
                    return;
                }
            }
            st = i__ + 1;
        }
        /* L60: */
    }
    /* Apply the singular values and treat the tiny ones as zero. */
    tol = rcnd * (r__1 = d__[aocl_blas_isamax(n, &d__[1], &c__1)], f2c_abs(r__1));
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Some of the elements in D can be negative because 1-by-1 */
        /* subproblems were not solved explicitly. */
        if((r__1 = d__[i__], f2c_abs(r__1)) <= tol)
        {
            aocl_lapack_slaset("A", &c__1, nrhs, &c_b6, &c_b6, &work[bx + i__ - 1], n);
        }
        else
        {
            ++(*rank);
            aocl_lapack_slascl("G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs,
                               &work[bx + i__ - 1], n, info);
        }
        d__[i__] = (r__1 = d__[i__], f2c_abs(r__1));
        /* L70: */
    }
    /* Now apply back the right singular vectors. */
    icmpq2 = 1;
    i__1 = nsub;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        st = iwork[i__];
        st1 = st - 1;
        nsize = iwork[sizei + i__ - 1];
        bxst = bx + st1;
        if(nsize == 1)
        {
            aocl_blas_scopy(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
        }
        else if(nsize <= *smlsiz)
        {
            aocl_blas_sgemm("T", "N", &nsize, nrhs, &nsize, &c_b11, &work[vt + st1], n, &work[bxst],
                            n, &c_b6, &b[st + b_dim1], ldb);
        }
        else
        {
            aocl_lapack_slalsa(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + b_dim1], ldb,
                               &work[u + st1], n, &work[vt + st1], &iwork[k + st1],
                               &work[difl + st1], &work[difr + st1], &work[z__ + st1],
                               &work[poles + st1], &iwork[givptr + st1], &iwork[givcol + st1], n,
                               &iwork[perm + st1], &work[givnum + st1], &work[c__ + st1],
                               &work[s + st1], &work[nwork], &iwork[iwk], info);
            if(*info != 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
        }
        /* L80: */
    }
    /* Unscale and sort the singular values. */
    aocl_lapack_slascl("G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, info);
    aocl_lapack_slasrt("D", n, &d__[1], info);
    aocl_lapack_slascl("G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], ldb, info);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLALSD */
}
/* slalsd_ */
