/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_hetrf_rook.c
 *  @brief Defines validate function of HETRF_ROOK() to use in test suite.
 *  */

#include "test_common.h"

void validate_hetrf_rook(char *uplo, integer n, integer lda, void *A_res, integer datatype,
                         integer *ipiv, double *residual, integer *info, void *A)
{
    if(n == 0)
    {
        return;
    }
    void *work = NULL;
    void *D = NULL;
    void *A_val = NULL;
    void *temp = NULL;
    void *X = NULL;
    void *B = NULL;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &D, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_val, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &work, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &temp, n);

    create_vector(datatype, &X, n);
    create_vector(datatype, &B, n);

    reset_matrix(datatype, n, n, D, n);

    set_identity_matrix(datatype, n, n, temp, n);
    set_identity_matrix(datatype, n, n, work, n);

    rand_vector(datatype, n, X, 1, d_zero, d_zero, 'R');

    integer s = 1, k;

    if(uplo[0] == 'U' || uplo[0] == 'u')
    {
        for(k = n - 1; k >= 0; k = k - s)
        {
            /*
             * Form U(K),
                      (   I    v    0   )   k-s
               U(k) = (   0    I    0   )   s
                      (   0    0    I   )   n-k
                         k-s   s   n-k
             */
            set_identity_matrix(datatype, n, n, A_val, n);
            /* U is a product of terms P(k) * U(k) where k decreases from n to 1 in steps of 1 or 2.
             * P(k) is a permutation matrix as defined by ipiv(k) */
            /* if ipiv(k) < 0 && ipiv(k-1) < 0, then rows k and -ipiv(k) of A_val are exchanged and
             * rows k-1 and -ipiv(k-1) of A_val are interchanged, D(k-1:k,k-1:k) is a 2-by-2
             * diagonal block */
            if(ipiv[k] < 0 && k > 0 && ipiv[k - 1] < 0)
            {
                /* Blocked Diagonal */
                s = 2;
                /* Get v matrix from API output */
                copy_submatrix(datatype, k - 1, 2, A_res, lda, A_val, n, 0, k - 1, 0, k - 1);
                /* Get blocked upper diagonal matrix from API output*/
                copy_submatrix(datatype, 2, 2, A_res, lda, D, n, k - 1, k - 1, k - 1, k - 1);
                copy_subvector(datatype, 1, A_res, lda, D, n, k - 1, k, k, k - 1);

                /* Diagonal blocks are hermitian. Set imaginary part of lower off diagonal element
                 * as negative of itself */
                negate_off_diagonal_element_imag(datatype, D, n, k, LOWER_OFF_DIAGONAL_ELEMENT);
                /* Swapping rows according to ipiv */
                swap_row_col(datatype, &n, A_val, n, &n, &n, (-1 * (ipiv[k - 1])) - 1, 0, k - 1, 0);
                swap_row_col(datatype, &n, A_val, n, &n, &n, (-1 * (ipiv[k])) - 1, 0, k, 0);
            }
            /* if ipiv(k) > 0, then rows and columns k and ipiv(k) of A_val are
             interchanged and D(k,k) is a 1-by-1 diagonal block */
            else
            {
                s = 1;
                /* Get v matrix from API output */
                copy_subvector(datatype, k, A_res, lda, A_val, n, 0, k, 0, k);
                /* Get diagonal matrix */
                copy_subvector(datatype, 1, A_res, lda, D, n, k, k, k, k);

                /* Swapping rows according to ipiv */
                if(ipiv[k] - 1 != k)
                {
                    swap_row_col(datatype, &n, A_val, n, &n, &n, ipiv[k] - 1, 0, k, 0);
                }
            }
            copy_matrix(datatype, "full", n, n, temp, n, work, n);
            /* Forming triangular matrix U from U(k) * U(k-1) * U(k-2) ... */
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, work, &n, A_val, &n, temp, &n);
        }
    }
    else if(uplo[0] == 'L' || uplo[0] == 'l')
    {
        for(k = 0; k <= n - 1; k += s)
        {
            /* Form L(k)
                       (   I    0     0   )  k-1
               L(k) =  (   0    I     0   )  s
                       (   0    v     I   )  n-k-s+1
                          k-1   s  n-k-s+1
            */
            set_identity_matrix(datatype, n, n, A_val, n);
            /* L is a product of terms P(k) * L(k) where k increases from 1 to n in steps of 1 or 2.
             * P(k) is a permutation matrix as defined by ipiv(k) */
            /* If ipiv(k) < 0 and ipiv(k+1) < 0, then rows k and -ipiv(k) of A_val are interchanged
             * and rows k+1 and -ipiv(k+1) of A_val are interchaged, D(k:k+1,k:k+1) is a 2-by-2
             * diagonal block */
            if(k < n - 1 && ipiv[k] < 0 && ipiv[k + 1] < 0)
            {
                /* Blocked Diagonal */
                s = 2;
                /* Get v matrix from API output */
                copy_submatrix(datatype, n - k - 2, 2, A_res, lda, A_val, n, k + s, k, k + s, k);
                /* Get blocked lower diagonal matrix from API output*/
                copy_submatrix(datatype, 2, 2, A_res, lda, D, n, k, k, k, k);
                copy_subvector(datatype, 1, A_res, lda, D, n, k + 1, k, k, k + 1);
                /* Diagonal blocks are hermitian. Set imaginary part of higher off diagonal element
                 * as negative of itself */
                negate_off_diagonal_element_imag(datatype, D, n, k, HIGHER_OFF_DIAGONAL_ELEMENT);
                /* Swapping rows according to ipiv */
                swap_row_col(datatype, &n, A_val, n, &n, &n, (-1 * (ipiv[k + 1])) - 1, 0, k + 1, 0);
                swap_row_col(datatype, &n, A_val, n, &n, &n, (-1 * (ipiv[k])) - 1, 0, k, 0);
            }
            /* If ipiv(k) > 0, then rows k and ipiv(k) of A_val are interchanged and D(k,k) is a
             * 1-by-1 diagonal block */
            else
            {
                s = 1;

                /* Get v matrix from API output */
                copy_subvector(datatype, n - k - 1, A_res, lda, A_val, n, k + 1, k, k + 1, k);
                /* Get diagonal matrix */
                copy_subvector(datatype, 1, A_res, lda, D, n, k, k, k, k);

                /* Swapping rows according to ipiv */
                if((ipiv[k] - 1) != k)
                {
                    swap_row_col(datatype, &n, A_val, n, &n, &n, (ipiv[k] - 1), 0, k, 0);
                }
            }
            copy_matrix(datatype, "full", n, n, temp, n, work, n);

            /* Forming triangular matrix L from L(1) * L(2) * L(3) ... */
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, work, &n, A_val, &n, temp, &n);
        }
    }
    switch(datatype)
    {
        case COMPLEX:
        {
            float norm_a, eps, resid1, resid2, norm;
            eps = fla_lapack_slamch("E");
            /* Test-1
             *Compute norm(A'*B - X)/(norm(X) * eps * n)
             */
            cgemv_("N", &n, &n, &c_one, A, &lda, X, &i_one, &c_zero, B, &i_one);
            chetrs_rook_(uplo, &n, &i_one, A_res, &lda, ipiv, B, &n, info);
            norm_a = fla_lapack_clange("1", &n, &i_one, X, &i_one, NULL);
            caxpy_(&n, &c_n_one, X, &i_one, B, &i_one);
            norm = fla_lapack_clange("1", &n, &i_one, B, &i_one, NULL);
            resid1 = norm / (eps * norm_a * n);

            /* Test-2
             * Compute norm(A-(U*D*U**T))/(norm(A) * eps * n)
             */
            if(*info > 0)
            {
                ((scomplex *)D)[(*info) * n + *info].real = 0;
                ((scomplex *)D)[(*info) * n + *info].imag = 0;
            }
            norm_a = fla_lapack_clange("1", &n, &n, A, &lda, NULL);
            cgemm_("N", "N", &n, &n, &n, &c_one, temp, &n, D, &n, &c_zero, A_val, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, A_val, &n, temp, &n, &c_n_one, A, &lda);
            norm = fla_lapack_clange("1", &n, &n, A, &lda, NULL);
            resid2 = norm / (eps * norm_a * n);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_a, eps, resid1, resid2, norm;
            eps = fla_lapack_dlamch("E");
            /* Test-1
             * Compute norm(A'*B - X)/(norm(X) * eps * n)
             */
            zgemv_("N", &n, &n, &z_one, A, &lda, X, &i_one, &z_zero, B, &i_one);
            zhetrs_rook_(uplo, &n, &i_one, A_res, &lda, ipiv, B, &n, info);
            norm_a = fla_lapack_zlange("1", &n, &i_one, X, &i_one, NULL);
            zaxpy_(&n, &z_n_one, X, &i_one, B, &i_one);
            norm = fla_lapack_zlange("1", &n, &i_one, B, &i_one, NULL);
            resid1 = norm / (eps * norm_a * n);

            /* Test-2
             * Compute norm(A-(U*D*U**T))/(norm(A) * eps * n)
             */
            if(*info > 0)
            {
                ((dcomplex *)D)[(*info) * n + *info].real = 0;
                ((dcomplex *)D)[(*info) * n + *info].imag = 0;
            }
            norm_a = fla_lapack_zlange("1", &n, &n, A, &lda, NULL);
            zgemm_("N", "N", &n, &n, &n, &z_one, temp, &n, D, &n, &z_zero, A_val, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, A_val, &n, temp, &n, &z_n_one, A, &lda);
            norm = fla_lapack_zlange("1", &n, &n, A, &lda, NULL);
            resid2 = norm / (eps * norm_a * n);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
    free_vector(B);
    free_vector(X);
    free_matrix(A_val);
    free_matrix(temp);
    free_matrix(D);
    free_matrix(work);
}
