/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gtsv.c
 *  @brief Defines validate function of GTSV() to use in test suite.
 *  */

#include "test_common.h"
void validate_gtsv(integer datatype, integer n, integer nrhs, void *B, integer ldb, void *X,
                   void *Xact, integer ldx, void *dl, void *d, void *du, void *dl_save,
                   void *d_save, void *du_save, integer info, void *scal, char imatrix, double *residual)
{
    if(n == 0 || nrhs == 0)
        return;

    void *work = NULL;
    void *A, *ipiv;
    integer dl_size, du_size, inc, lda;

    /* Set sizes and dimensions for vector and matrix */
    dl_size = n - 2, du_size = n - 1;
    lda = ldb;
    inc = lda + 1;

    create_vector(INTEGER, &ipiv, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, ldb);
    /* A matrix consist of the tridiagonal input vectors */
    copy_tridiag_matrix(datatype, dl, d, du, n, n, A, ldb);
    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm_X = 0, norm_B = 0, norm1 = 0, norm2 = 0, resid1, resid2;
            float norm_dl = 0, norm_d = 0, norm_du = 0, norm3 = 0, norm4 = 0, norm5 = 0, resid3;
            float *d_ptr, *du_ptr, *dl_ptr;

            /* Test1: Compute norm(AX-B / (eps * norm_B * n)) */
            eps = fla_lapack_slamch("P");
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    float *vector = (float *)B + i * ldb;
                    norm_B = fla_max(norm_B, snrm2_(&n, vector, &i_one));
                }
            }
            else
                norm_B = fla_lapack_slange("F", &n, &nrhs, B, &ldb, work);
            slagtm_("N", &n, &nrhs, &s_one, dl, d, du, X, &ldx, &s_n_one, B, &ldb);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    float *vector = (float *)B + i * ldb;
                    norm1 = fla_max(norm1, snrm2_(&n, vector, &i_one));
                }
            }
            else
                norm1 = fla_lapack_slange("F", &n, &nrhs, B, &ldb, work);
            resid1 = (norm1 / norm_B)/(eps * n);

            /* Test 2: Compute norm (Xact - X / (eps * norm_X *n)) */
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    float *vector = (float *)X + i * ldx;
                    norm_X = fla_max(norm_X, snrm2_(&n, vector, &i_one));
                }
            }
            else
                norm_X = fla_lapack_slange("F", &n, &nrhs, X, &ldx, work);
            matrix_difference(datatype, n, nrhs, X, ldx, Xact, ldx);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    float *vector = (float *)X + i * ldx;
                    norm2 = fla_max(norm2, snrm2_(&n, vector, &i_one));
                }
            }
            else
                norm2 = fla_lapack_slange("F", &n, &nrhs, X, &ldx, work);
            resid2 = (norm2 / norm_X) /(eps * n);

            /* Test3 : Compute norm for tridiagonal vectors dl, d, du by
             * Calculating the difference between tridiagonal vector (d, du, dl)
             * and upper triangular matrix of A from getrf */
            norm_dl = fla_lapack_slange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm_d = fla_lapack_slange("F", &n, &i_one, d_save, &i_one, work);
            norm_du = fla_lapack_slange("F", &du_size, &i_one, du_save, &i_one, work);
            sgetrf_(&n, &n, A, &lda, ipiv, &info);
            /* Calculate the difference between d_save and the diagonal of A matrix */
            d_ptr = ((float *)A);
            saxpy_(&n, &s_n_one, d_ptr, &inc, d_save, &i_one);
            /* Calculate the difference between du_save and the upper diagonal of A matrix*/
            du_ptr = ((float *)A + lda);
            saxpy_(&du_size, &s_n_one, du_ptr, &inc, du_save, &i_one);
            /* Calculate the difference between dl_save and the 2nd upper diagonal of A matrix */
            dl_ptr = ((float *)A + (2 * lda));
            saxpy_(&dl_size, &s_n_one, dl_ptr, &inc, dl_save, &i_one);

            norm3 = fla_lapack_slange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm4 = fla_lapack_slange("F", &n, &i_one, d_save, &i_one, work);
            norm5 = fla_lapack_slange("F", &du_size, &i_one, du_save, &i_one, work);
            resid3 = ((norm3 + norm4 + norm5) / (norm_dl + norm_d + norm_du)) / (eps * n);

            *residual = (double)fla_max(fla_max(resid1, resid2), resid3);
            break;
        }
        case DOUBLE:
        {
            double eps, norm_X = 0, norm_B = 0, norm1 = 0, norm2 = 0, resid1 = 0, resid2 = 0;
            double norm_dl = 0, norm_d = 0, norm_du = 0, norm3 = 0, norm4 = 0, norm5 = 0, resid3 = 0;
            double *d_ptr, *dl_ptr, *du_ptr;

            /* Test1: Compute norm(AX-B / (eps * norm(b) * n)) */
            eps = fla_lapack_dlamch("P");
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    double *vector = (double *)B + i * ldb;
                    norm_B = fla_max(norm_B, dnrm2_(&n, vector, &i_one));
                }
            }
            else
                norm_B = fla_lapack_dlange("F", &n, &nrhs, B, &ldb, work);
            dlagtm_("N", &n, &nrhs, &d_one, dl, d, du, X, &ldx, &d_n_one, B, &ldb);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    double *vector = (double *)B + i * ldb;
                    norm1 = fla_max(norm1, dnrm2_(&n, vector, &i_one));
                }
            }
            else
                norm1 = fla_lapack_dlange("F", &n, &nrhs, B, &ldb, work);
            resid1 = (norm1 / norm_B) /(eps * n);

            /* Test 2: Compute norm (xact - x / (eps * norm(x) *n)) */
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    double *vector = (double *)X + i * ldx;
                    norm_X = fla_max(norm_X, dnrm2_(&n, vector, &i_one));
                }
            }
            else
                norm_X = fla_lapack_dlange("F", &n, &nrhs, X, &ldx, work);
            matrix_difference(datatype, n, nrhs, X, ldx, Xact, ldx);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    double *vector = (double *)X + i * ldx;
                    norm2 = fla_max(norm2, dnrm2_(&n, vector, &i_one));
                }
            }
            else
                norm2 = fla_lapack_dlange("F", &n, &nrhs, X, &ldx, work);
            resid2 = (norm2 / norm_X) /(eps * n);

            /* Test3 : Compute norm for tridiagonal vectors dl, d, du by
             * Calculating the difference between tridiagonal vector (d, du, dl)
             * and upper triangular matrix of A from getrf */
            norm_dl = fla_lapack_dlange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm_d = fla_lapack_dlange("F", &n, &i_one, d_save, &i_one, work);
            norm_du = fla_lapack_dlange("F", &du_size, &i_one, du_save, &i_one, work);
            dgetrf_(&n, &n, A, &lda, ipiv, &info);
            /* Calculate the difference between d_save and the diagonal of A matrix */
            d_ptr = ((double *)A);
            daxpy_(&n, &d_n_one, d_ptr, &inc, d_save, &i_one);
            /* Calculate the difference between du_save and the upper diagonal of A matrix */
            du_ptr = ((double *)A + lda);
            daxpy_(&du_size, &d_n_one, du_ptr, &inc, du_save, &i_one);
            /* Calculate the difference between dl_save and the 2nd upper diagonal of A matrix */
            dl_ptr = ((double *)A + (2 * lda));
            daxpy_(&dl_size, &d_n_one, dl_ptr, &inc, dl_save, &i_one);

            norm3 = fla_lapack_dlange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm4 = fla_lapack_dlange("F", &n, &i_one, d_save, &i_one, work);
            norm5 = fla_lapack_dlange("F", &du_size, &i_one, du_save, &i_one, work);
            resid3 = ((norm3 + norm4 + norm5) / (norm_dl + norm_d + norm_du)) / (eps * n);

            *residual = (double)fla_max(fla_max(resid1, resid2), resid3);
            break;
        }
        case COMPLEX:
        {
            float eps, norm_X = 0, norm_B = 0, norm1 = 0, norm2 = 0, resid1 = 0, resid2;
            float norm_dl = 0, norm_d = 0, norm_du = 0, norm3 = 0, norm4 = 0, norm5 = 0, resid3;
            scomplex *d_ptr, *dl_ptr, *du_ptr;
            /* Test1: Compute norm(AX-B / (eps * norm_B * n)) */
            eps = fla_lapack_slamch("P");
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    scomplex *vector = (scomplex *)B + i * ldb;
                    norm_B = fla_max(norm_B, scnrm2_(&n, vector, &i_one));
                }
            }
            else
                norm_B = fla_lapack_clange("F", &n, &nrhs, B, &ldb, work);
            clagtm_("N", &n, &nrhs, &s_one, dl, d, du, X, &ldx, &s_n_one, B, &ldb);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    scomplex *vector = (scomplex *)B + i * ldb;
                    norm1 = fla_max(norm1, scnrm2_(&n, vector, &i_one));
                }
            }
            else
                norm1 = fla_lapack_clange("F", &n, &nrhs, B, &ldb, work);
            resid1 = (norm1 / norm_B) /(eps * n);

            /* Test 2: Compute norm (xact - x / (eps * norm_X *n)) */
            norm_X = fla_lapack_clange("F", &n, &nrhs, X, &ldx, work);
            matrix_difference(datatype, n, nrhs, X, ldx, Xact, ldx);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    scomplex *vector = (scomplex *)X + i * ldx;
                    norm_X = fla_max(norm_X, scnrm2_(&n, vector, &i_one));
                }
            }
            else
                norm2 = fla_lapack_clange("F", &n, &nrhs, X, &ldx, work);
            resid2 = (norm2 / norm_X) /(eps * n);

            /* Test3 : Compute norm for tridiagonal vectors dl, d, du by
             * Calculating the difference between tridiagonal vector (d, du, dl)
             * and upper triangular matrix of A from getrf */
            norm_dl = fla_lapack_clange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm_d = fla_lapack_clange("F", &n, &i_one, d_save, &i_one, work);
            norm_du = fla_lapack_clange("F", &du_size, &i_one, du_save, &i_one, work);
            cgetrf_(&n, &n, A, &lda, ipiv, &info);
            /* Calculate the difference between d_save and the diagonal of A matrix */
            d_ptr = ((scomplex *)A);
            caxpy_(&n, &c_n_one, d_ptr, &inc, d_save, &i_one);
            /* Calculate the difference between du_save and the upper diagonal of A matrix */
            du_ptr = ((scomplex *)A + lda);
            caxpy_(&du_size, &c_n_one, du_ptr, &inc, du_save, &i_one);
            /* Calculate the difference between dl_save and the 2nd upper diagonal of A matrix */
            dl_ptr = ((scomplex *)A + (2 * lda));
            caxpy_(&dl_size, &c_n_one, dl_ptr, &inc, dl_save, &i_one);

            norm3 = fla_lapack_clange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm4 = fla_lapack_clange("F", &n, &i_one, d_save, &i_one, work);
            norm5 = fla_lapack_clange("F", &du_size, &i_one, du_save, &i_one, work);
            resid3 = ((norm3 + norm4 + norm5) / (norm_dl + norm_d + norm_du)) / (eps * n);

            *residual = (double)fla_max(fla_max(resid1, resid2), resid3);

            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_X = 0, norm_B = 0, norm1 = 0, norm2 = 0, eps, resid1, resid2;
            double norm_dl = 0, norm_d = 0, norm_du = 0, norm3 = 0, norm4 = 0, norm5 = 0, resid3;
            dcomplex *d_ptr, *dl_ptr, *du_ptr;

            /* Test1: Compute norm(AX-B / (eps * norm_B * n)) */
            eps = fla_lapack_dlamch("P");
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    dcomplex *vector = (dcomplex *)B + i * ldb;
                    norm_B = fla_max(norm_B, dznrm2_(&n, vector, &i_one));
                }
            }
            else
                norm_B = fla_lapack_zlange("F", &n, &nrhs, B, &ldb, work);
            zlagtm_("N", &n, &nrhs, &d_one, dl, d, du, X, &ldx, &d_n_one, B, &ldb);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    dcomplex *vector = (dcomplex *)B + i * ldb;
                    norm1 = fla_max(norm1, dznrm2_(&n, vector, &i_one));
                }
            }
            else
                norm1 = fla_lapack_zlange("F", &n, &nrhs, B, &ldb, work);
            resid1 = (norm1 / norm_B) /(eps * n);

            /* Test 2: Compute norm (xact - x / (eps * norm_X *n)) */
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    dcomplex *vector = (dcomplex *)X + i * ldx;
                    norm_X = fla_max(norm_X, dznrm2_(&n, vector, &i_one));
                }
            }
            else
                norm_X = fla_lapack_zlange("F", &n, &nrhs, X, &ldx, work);
            matrix_difference(datatype, n, nrhs, X, ldx, Xact, ldx);
            if(imatrix == 'O')
            {
                for(int i = 0; i < nrhs; i++)
                {
                    /* To handle large size values nrm2 is used */
                    dcomplex *vector = (dcomplex *)X + i * ldx;
                    norm2 = fla_max(norm2, dznrm2_(&n, vector, &i_one));
                }
            }
            else
                norm2 = fla_lapack_zlange("F", &n, &nrhs, X, &ldx, work);
            resid2 = (norm2 / norm_X) /(eps * n);

            /* Test3 : Compute norm for tridiagonal vectors dl, d, du by
             * Calculating the difference between tridiagonal vector (d, du, dl)
             * and upper triangular matrix of A from getrf */
            norm_dl = fla_lapack_zlange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm_d = fla_lapack_zlange("F", &n, &i_one, d_save, &i_one, work);
            norm_du = fla_lapack_zlange("F", &du_size, &i_one, du_save, &i_one, work);
            zgetrf_(&n, &n, A, &lda, ipiv, &info);
            /* Calculate the difference between d_save and the diagonal of A matrix */
            d_ptr = ((dcomplex *)A);
            zaxpy_(&n, &z_n_one, d_ptr, &inc, d_save, &i_one);
            /* Calculate the difference between du_save and the upper diagonal of A matrix */
            du_ptr = ((dcomplex *)A + lda);
            zaxpy_(&du_size, &z_n_one, du_ptr, &inc, du_save, &i_one);
            /* Calculate the difference between dl_save and the 2nd upper diagonal of A matrix */
            dl_ptr = ((dcomplex *)A + (2 * lda));
            zaxpy_(&dl_size, &z_n_one, dl_ptr, &inc, dl_save, &i_one);

            norm3 = fla_lapack_zlange("F", &dl_size, &i_one, dl_save, &i_one, work);
            norm4 = fla_lapack_zlange("F", &n, &i_one, d_save, &i_one, work);
            norm5 = fla_lapack_zlange("F", &du_size, &i_one, du_save, &i_one, work);
            resid3 = ((norm3 + norm4 + norm5) / (norm_dl + norm_d + norm_du)) / (eps * n);

            *residual = (double)fla_max(fla_max(resid1, resid2), resid3);
            break;
        }
    }
    free_vector(ipiv);
    free_matrix(A);
}
