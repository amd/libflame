/******************************************************************************
* Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gels.c
 *  @brief Defines validate function of GELS() to use in test suite.
 *  */

#include "test_common.h"

void validate_gels(char *trans,
                   integer m,
                   integer n,
                   integer nrhs,
                   void *A,
                   integer lda,
                   void *B,
                   integer ldb,
                   void *x,
                   integer datatype,
                   double *residual,
                   integer *info)
{
    integer m1 = m, n1 = n, INFO = 0, lwork = (m+nrhs)*(n+2), ldwork = (m + nrhs), temp;
    char NORM = '1';
    void *work = NULL;
    void *C = NULL;
    double temp1;

    if(*trans == 'T' || *trans == 'C')
    {
        m1 = n;
        n1 = m;
        NORM = 'I';
        lwork = (n+nrhs)*(m+2);
        ldwork = m;
    }

    create_matrix(datatype, &C, ldb, n1);
    create_vector(datatype, &work, lwork);

    if(m == 0 || n == 0)
    {
        return;
    }
    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            float resid1 = 0, resid2 = 0, resid3 = 0, rwork;
            eps = fla_lapack_slamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'T' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                resid1 = eps;
                temp = m1-n1;
                for(int i = 0; i < nrhs; i++)
                {
                    resid1 = snrm2_(&temp, &((float*)x)[(i * ldb) + n1], &i_one);
                    resid1 = fla_max(norm, resid1);
                }
                *residual = (double)resid1;
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_slange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_slange(&NORM, &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_slange(&NORM, &n1, &nrhs, x, &i_one, work);
                sgemm_(trans, "N", &m1, &nrhs, &n1, &s_n_one, A, &lda, x, &ldb, &s_one, B, &ldb);
                norm = fla_lapack_slange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1 ,n1 ) * norm_a * norm_x * eps);

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                sgemm_("T", "N", &nrhs, &n1, &m1, &s_one, B, &ldb, A, &lda, &s_zero, C, &ldb);
                norm = fla_lapack_slange("M", &nrhs, &n1, C, &ldb, work);
                resid2 = norm / (fla_max(m1 ,fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                /* Test - 3
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 */
                if((*trans == 'T' && m > n) || (*trans == 'N' && m < n))
                {
                    /* Copy A into work */
                    fla_lapack_slacpy("All", &m, &n, A, &lda, work, &ldwork);
                    norm_a = fla_lapack_slange("M", &m, &n, work, &ldwork, &rwork);
                    /*Scale work*/
                    if (norm_a != 0.)
                    {
                        slascl_("G", &i_zero, &i_zero, &norm_a, &s_one, &m, &n, work, &ldwork,
                                &INFO);
                    }
                    if(*trans == 'T')
                    {
                        /*Copy x into work*/
                        fla_lapack_slacpy("All", &m, &nrhs, x, &ldb,
                                          &((float*)work)[n * ldwork], &ldwork);
                        norm_x = fla_lapack_slange("M", &m, &nrhs, &((float*)work)[n * ldwork],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                            slascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &m, &nrhs,
                                    &((float*)work)[n * ldwork], &ldwork, &INFO);
                        }
                        temp = n + nrhs;
                        /*QR factorization of x*/
                        sgeqr2_(&m, &temp, work, &ldwork, &((float*)work)[ldwork * (n + nrhs)],
                                &((float*)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))],
                                &INFO);
                        norm = 0;
                        /*Compute norm*/
                        for(int j = n; j <= temp; j++)
                        {
                            for(int i = n; i < fla_min(m, j); i++)
                            {
                                temp1 = FLA_FABS(((float*)work)[i + (j - 1) * m]);
                                norm = fla_max(temp1, norm);
                            }
                        }
                    }
                    else if( *trans == 'N')
                    {
                        /*Copy x into work*/
                        for( int i = 0; i < n; i++)
                        {
                            for(int j = 0; j < nrhs; j++)
                            {
                                ((float*)work)[m + j + (i * ldwork)] = ((float*)x)[i + j * ldb];
                            }
                        }
                        norm_x = fla_lapack_slange("M", &nrhs, &n, &((float*)work)[m],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                        	slascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &nrhs, &n,
                                    &((float*)work)[m + 1], &ldwork, &INFO);
                        }
                        /*LQ factorization*/
                        sgelq2_(&ldwork, &n, work, &ldwork, &((float*)work)[ldwork * n],
                                &((float*)work)[ldwork * (n + 1)], &INFO);
                        /*Compute norm*/
                        for(int j = n; j <= n; j++)
                        {
                            for(int i = n; i < ldwork; i++)
                            {
                                temp1 = FLA_FABS(((float*)work)[i + (j - 1) * ldwork]);
                                norm = fla_max(temp1, norm);
                            }
                        }
                    }
                    resid3 = norm / ((double) fla_max(m, fla_max(n, nrhs)) * eps);
                }
                *residual = (double)fla_max(resid1, fla_max(resid2, resid3));
            }
            break;
        }
        case DOUBLE:
        {
            double eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            double resid1 = 0, resid2 = 0, resid3 = 0, rwork;
            eps = fla_lapack_dlamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'T' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                resid1 = eps;
                temp = m1 - n1;
                for(int i = 0; i < nrhs; i++)
                {
                    resid1 = dnrm2_(&temp, &((double*)x)[(i * ldb) + n1], &i_one);
                    resid1 = fla_max(norm, resid1);
                }
                *residual = (double)resid1;
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_dlange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_dlange(&NORM, &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_dlange(&NORM, &n1, &nrhs, x, &i_one, work);
                dgemm_(trans, "N", &m1, &nrhs, &n1, &d_n_one, A, &lda, x, &ldb, &d_one, B, &ldb);
                norm = fla_lapack_dlange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1 ,n1 ) * norm_a * norm_x * eps);

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                dgemm_("T", "N", &nrhs, &n1, &m1, &d_one, B, &ldb, A, &lda, &d_zero, C, &ldb);
                norm = fla_lapack_dlange("1", &nrhs, &n1, C, &ldb, work);
                resid2 = norm / (fla_max(m1 ,fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                /* Test - 3
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 */
                if((*trans == 'T' && m > n) || (*trans == 'N' && m < n))
                {
                    /*Copy A into work*/
                    fla_lapack_dlacpy("All", &m, &n, A, &lda, work, &ldwork);
                    norm_a = fla_lapack_dlange("M", &m, &n, work, &ldwork, &rwork);
                    /*Scale work*/
                    if (norm_a != 0.)
                    {
                        dlascl_("G", &i_zero, &i_zero, &norm_a, &d_one, &m, &n, work, &ldwork,
                                &INFO);
                    }
                    if(*trans == 'T')
                    {
                        /*Copy x into work*/
                        fla_lapack_dlacpy("All", &m, &nrhs, x, &ldb,
                                          &((double*)work)[n * ldwork], &ldwork);
                        norm_x = fla_lapack_dlange("M", &m, &nrhs, &((double*)work)[n * ldwork],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                            dlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &m, &nrhs,
                                    &((double*)work)[n * ldwork], &ldwork, &INFO);
                        }
                        temp = n + nrhs;
                        /*QR factorization of x*/
                        dgeqr2_(&m, &temp, work, &ldwork, &((double*)work)[ldwork * (n + nrhs)],
                                &((double*)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))],
                                &INFO);
                        norm = 0;
                        /*Compute norm*/
                        for(int j = n; j <= temp; j++)
                        {
                            for(int i = n; i < fla_min(m, j); i++)
                            {
                                temp1 = FLA_FABS(((double*)work)[i + (j - 1) * m]);
                                norm = fla_max(temp1, norm);
                            }
                        }
                    }
                    else if( *trans == 'N')
                    {
                        /*Copy x into work*/
                        for( int i = 0; i < n; i++)
                        {
                            for(int j = 0; j < nrhs; j++)
                            {
                                ((double*)work)[m + j + (i * ldwork)] = ((double*)x)[i + j * ldb];
                            }
                        }
                        norm_x = fla_lapack_dlange("M", &nrhs, &n, &((double*)work)[m],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                            dlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &nrhs, &n,
                                    &((double*)work)[m + 1], &ldwork, &INFO);
                        }
                        /*LQ factorization*/
                        dgelq2_(&ldwork, &n, work, &ldwork, &((double*)work)[ldwork * n],
                                &((double*)work)[ldwork * (n + 1)], &INFO);
                        /*Compute norm*/
                        for(int j = n; j <= n; j++)
                        {
                            for(int i = n; i < ldwork; i++)
                            {
                                temp1 = FLA_FABS(((double*)work)[i + (j - 1) * ldwork]);
                                norm = fla_max(temp1, norm);
                            }
                        }
                    }
                    resid3 = norm / ((double) fla_max(m, fla_max(n, nrhs)) * eps);
                }
                *residual = (double)fla_max(resid1, fla_max(resid2, resid3));
            }
            break;
        }
        case COMPLEX:
        {
            float eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            float resid1 = 0, resid2 = 0, resid3 = 0, rwork;
            eps = fla_lapack_slamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'C' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                resid1 = eps;
                temp = m1 - n1;
                for(int i = 0; i < nrhs; i++)
                {
                    resid1 = scnrm2_(&temp, &((scomplex*)x)[(i * ldb) + n1], &i_one);
                    resid1 = fla_max(norm, resid1);
                }
                *residual = (double)resid1;
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_clange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_clange(&NORM, &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_clange(&NORM, &n1, &nrhs, x, &i_one, work);
                cgemm_(trans, "N", &m1, &nrhs, &n1, &c_n_one, A, &lda, x, &ldb, &c_one, B, &ldb);
                norm = fla_lapack_clange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1 ,n1 ) * norm_a * norm_x * eps);
                *residual = (double)resid1;

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                cgemm_("T", "N", &nrhs, &n1, &m1, &c_one, B, &ldb, A, &lda, &c_zero, C, &ldb);
                norm = fla_lapack_clange("1", &nrhs, &n1, C, &ldb, work);
                resid2 = norm / (fla_max(m1 ,fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                /* Test - 3
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 */
                if((*trans == 'C' && m > n) || (*trans == 'N' && m < n))
                {
                    /*Copy A into work*/
                    fla_lapack_clacpy("All", &m, &n, A, &lda, work, &ldwork);
                    norm_a = fla_lapack_clange("M", &m, &n, work, &ldwork, &rwork);
                    /*Scale work*/
                    if (norm_a != 0.)
                    {
	                    clascl_("G", &i_zero, &i_zero, &norm_a, &s_one, &m, &n, work, &ldwork,
	                	       &INFO);
                    }
                    if(*trans == 'C')
                    {
                        /*Copy x into work*/
                        fla_lapack_clacpy("All", &m, &nrhs, x, &ldb,
                                          &((scomplex*)work)[n * ldwork], &ldwork);
                        norm_x = fla_lapack_clange("M", &m, &nrhs, &((scomplex*)work)[n * ldwork],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                            clascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &m, &nrhs,
                                    &((scomplex*)work)[n * ldwork], &ldwork, &INFO);
                        }
                        temp = n + nrhs;
                        /*QR factorization of x*/
                        cgeqr2_(&m, &temp, work, &ldwork, &((scomplex*)work)[ldwork * (n + nrhs)],
                                &((scomplex*)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))],
                                &INFO);
                        norm = 0;
                        /*Compute norm*/
                        for(int j = n; j <= temp; j++)
                        {
                            for(int i = n; i < fla_min(m, j); i++)
                            {
                                temp1 = FLA_FABS(((scomplex*)work)[i + (j - 1) * m].real);
                                norm = fla_max(temp1, norm);
                            }
                        }
                    }
                    else if( *trans == 'N')
                    {
                        /*Copy x into work*/
                        for( int i = 0; i < n; i++)
                        {
                            for(int j = 0; j < nrhs; j++)
                            {
                                ((scomplex*)work)[m + j + (i * ldwork)] = ((scomplex*)x)[i + j * ldb];
                            }
                        }
                        norm_x = fla_lapack_clange("M", &nrhs, &n, &((scomplex*)work)[m],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                            clascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &nrhs, &n,
                                    &((scomplex*)work)[m + 1], &ldwork, &INFO);
                        }
                        /*LQ factorization*/
                        cgelq2_(&ldwork, &n, work, &ldwork, &((scomplex*)work)[ldwork * n],
                                &((scomplex*)work)[ldwork * (n + 1)], &INFO);
                        /*Compute norm*/
                        for(int j = n; j <= n; j++)
                        {
                            for(int i = n; i < ldwork; i++)
                            {
                                temp1 = FLA_FABS(((scomplex*)work)[i + (j - 1) * ldwork].real);
                                norm = fla_max(norm, temp1);
                            }
                        }
                    }
                    resid3 = norm / ((double) fla_max(m, fla_max(n, nrhs)) * eps);

                }
                *residual = (double)fla_max(resid1, fla_max(resid2, resid3));
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            double resid1 = 0, resid2 = 0, resid3 = 0, rwork;
            eps = fla_lapack_slamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'C' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                resid1 = eps;
                temp = m1 - n1;
                for(int i = 0; i < nrhs; i++)
                {
                    resid1 = dznrm2_(&temp, &((dcomplex*)x)[(i * ldb) + n1], &i_one);
                    resid1 = fla_max(norm, resid1);
                }
                *residual = (double)resid1;
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_zlange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_zlange(&NORM, &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_zlange(&NORM, &n1, &nrhs, x, &i_one, work);
                zgemm_(trans, "N", &m1, &nrhs, &n1, &z_n_one, A, &lda, x, &ldb, &z_one, B, &ldb);
                norm = fla_lapack_zlange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1 ,n1 ) * norm_a * norm_x * eps);
                *residual = (double)resid1;

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                zgemm_("T", "N", &nrhs, &n1, &m1, &z_one, B, &ldb, A, &lda, &z_zero, C, &ldb);
                norm = fla_lapack_zlange("1", &nrhs, &n1, C, &ldb, work);
                resid2 = norm / (fla_max(m1 ,fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                /* Test - 3
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 */
                if((*trans == 'C' && m > n) || (*trans == 'N' && m < n))
                {
                    /*Copy A into work*/
                    fla_lapack_zlacpy("All", &m, &n, A, &lda, work, &ldwork);
                    norm_a = fla_lapack_zlange("M", &m, &n, work, &ldwork, &rwork);
                    /*Scale work*/
                    if (norm_a != 0)
                    {
                        zlascl_("G", &i_zero, &i_zero, &norm_a, &d_one, &m, &n, work, &ldwork,
                                &INFO);
                    }
                    if(*trans == 'C')
                    {
                        /*Copy x into work*/
                        fla_lapack_zlacpy("All", &m, &nrhs, x, &ldb,
                                          &((dcomplex*)work)[n * ldwork], &ldwork);
                        norm_x = fla_lapack_zlange("M", &m, &nrhs, &((dcomplex*)work)[n * ldwork],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                            zlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &m, &nrhs,
                                    &((dcomplex*)work)[n * ldwork], &ldwork, &INFO);
                        }
                        temp = n + nrhs;
                        /*QR factorization of x*/
                        zgeqr2_(&m, &temp, work, &ldwork, &((dcomplex*)work)[ldwork * (n + nrhs)],
                                &((dcomplex*)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))],
                                &INFO);
                        norm = 0;
                        /*Compute norm*/
                        for(int j = n; j <= temp; j++)
                        {
                            for(int i = n; i < fla_min(m, j); i++)
                            {
                                temp1 = FLA_FABS(((dcomplex*)work)[i + (j - 1) * m].real);
                                norm = fla_max(temp1, norm);
                            }
                        }
                    }
                    else if( *trans == 'N')
                    {
                        /*Copy x into work*/
                        for( int i = 0; i < n; i++)
                        {
                            for(int j = 0; j < nrhs; j++)
                            {
                                ((dcomplex*)work)[m + j + (i * ldwork)] = ((dcomplex*)x)[i + j * ldb];
                            }
                        }
                        norm_x = fla_lapack_zlange("M", &nrhs, &n, &((dcomplex*)work)[m],
                                                   &ldwork, &rwork);
                        /*Scale x*/
                        if (norm_x != 0)
                        {
                        	zlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &nrhs, &n,
                                    &((dcomplex*)work)[m + 1], &ldwork, &INFO);
                        }
                        /*LQ factorization*/
                        zgelq2_(&ldwork, &n, work, &ldwork, &((dcomplex*)work)[ldwork * n],
                                &((dcomplex*)work)[ldwork * (n + 1)], &INFO);
                        /*Compute norm*/
                        for(int j = n; j <= n; j++)
                        {
                            for(int i = n; i < ldwork; i++)
                            {
                                temp1 = FLA_FABS(((dcomplex*)work)[i + (j - 1) * ldwork].real);
                                norm = fla_max(norm, temp1);
                            }
                        }
                    }
                    resid3 = norm / ((double) fla_max(m, fla_max(n, nrhs)) * eps);
                }
                *residual = (double)fla_max(resid1, fla_max(resid2, resid3));

            }
            break;
        }
    }
    free_vector(work);
    free_matrix(C);
}
