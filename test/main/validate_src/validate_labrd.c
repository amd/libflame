/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_labrd.c
 *  @brief Defines validate function of LABRD() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;
integer Ut = 0, Vt = 0;

/* Helper functions for collecting output from LABRD and recreating A in bidiagonal form */
void collect_LABRD_Output(integer datatype, integer m, integer n, integer nb, void *A_test,
                          integer lda, void *X, integer ldx, void *Y, integer ldy, void *U, void *V,
                          void *XU, void *VY);
void create_bidiagonal_form(integer datatype, integer m, integer n, integer nb, void *A,
                            integer lda, void *XU, void *YV, void *d, void *e);

void validate_labrd(char *tst_api, integer m, integer n, integer nb, void *A, void *A_test,
                    integer lda, void *d, void *e, void *tauq, void *taup, void *X, integer ldx,
                    void *Y, integer ldy, integer datatype, double err_thresh, FILE *g_ext_fptr,
                    char imatrix)
{
    double residual = err_thresh;
    /* Early return conditions */
    if(m == 0 || n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }
    /* Print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);
    /*
     * For validation of LABRD API,
     * We first create the bidiagonal matrix using the outputs from LABRD.
     * We then create a secondary bidiagonal matrix by applying the elementary
     * reflectors in a sequential order following the theory of bidiagonal
     * reduction from
     * https://www.cs.utexas.edu/~flame/laff/alaff/chapter11-reduction-to-bidirectional-form.html
     * We then get the difference from output generated using the two methods to
     * get the error in output.
     *
     */
    void *V = NULL, *U = NULL, *XU = NULL, *VY = NULL;
    void *A_BD = NULL;
    /* Creating matrices for outputs to be collected from LABRD */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, nb, &V, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, nb, n, &U, nb);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &XU, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &VY, m);
    /* Creating a matrix which will be used to hold A in bidiagonal form
     * generated from LABRD outputs*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_BD, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_BD, lda);
    /* Vt and Ut are used for identifying if the subdiagonal is present in
     * the upper triangle or lower triangle of the reduced matrix */
    if(m >= n)
    {
        Vt = 0;
        Ut = 1;
    }
    else
    {
        Vt = 1;
        Ut = 0;
    }
    /*
     * Collect U and V matrices from output of labrd
     * If m >= n,
     * v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in A(i:m,i);
     * u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n)
     *
     * If m < n,
     * v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in A(i+2:m,i);
     * u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in A(i,i+1:n)
     *
     * The elements of the vectors v and u together form the m-by-nb matrix
     * V and the nb-by-n matrix U**H, Which are used to form X * U**H and
     * V * Y**H which are used in defining the bidiagonal matrix output from
     * LABRD
     *
     */
    collect_LABRD_Output(datatype, m, n, nb, A_test, lda, X, ldx, Y, ldy, U, V, XU, VY);
    /*
     * Creating bidiagonal reduced A from outputs of LABRD
     *
     *  - The first NB diagonal elements are defined by D output of LABRD
     *  - The first NB sub diagonal elements are defined by E output of LABRD
     *  - The sub diagonal being in the upper or lower triangle is decided
     * by the relation between m and n.
     *  - The values present in the first NB rows and columns apart from the
     * diagonals are zero.
     *  - The values present in the remaining positions for the bidiagonal matrix
     * are defined by (A - V*Y**T - X*U**T) which represents the reduced form of
     * the matrix
     *
     */
    create_bidiagonal_form(datatype, m, n, nb, A_BD, lda, XU, VY, d, e);
    /*
     * Create bidiagonal form by applying elementary reflectors
     *
     * The bidiagonal form of a matrix is obtained by sequentially applying
     * the column elementary reflectors to the left and row elementary reflectors
     * to the right of the matrix.
     * We generate another matrix with the first NB columns in bidiagonal form and
     * the remaining values in the reduced form by applying this method.
     *
     */
    switch(datatype)
    {
        case FLOAT:
        {
            for(integer i = 0; i < nb; i++)
            {
                void *LR = NULL, *RR = NULL, *work = NULL;
                float tauQ, tauP;
                /* Collecting elementary reflector V from outputs in A */
                tauQ = ((float *)tauq)[i];
                create_vector(datatype, &LR, m);
                for(int j = 0; j < m; j++)
                {
                    if(j < i + Vt)
                        ((float *)LR)[j] = s_zero;
                    else if(j == i + Vt)
                        ((float *)LR)[j] = s_one;
                    else
                        ((float *)LR)[j] = ((float *)A_test)[j + i * lda];
                }
                /* Applying elementary reflector V to the left of the matrix */ 
                create_vector(datatype, &work, n);
                slarf_("L", &m, &n, LR, &i_one, &tauQ, A, &lda, work);
                free_vector(LR);
                free_vector(work);
                /* Collecting elementary reflector U from outputs in A */
                tauP = ((float *)taup)[i];
                create_vector(datatype, &RR, n);
                for(int j = 0; j < n; j++)
                {
                    if(j < i + Ut)
                        ((float *)RR)[j] = s_zero;
                    else if(j == i + Ut)
                        ((float *)RR)[j] = s_one;
                    else
                        ((float *)RR)[j] = ((float *)A_test)[i + j * lda];
                }
                /* Applying elementary reflector U to the right of the matrix */
                create_vector(datatype, &work, m);
                slarf_("R", &m, &n, RR, &i_one, &tauP, A, &lda, work);
                free_vector(RR);
                free_vector(work);
            }
            /* Get difference between the 2 bidiagonal forms of A */
            float norm, norm_A, eps;
            eps = fla_lapack_slamch("P");
            norm_A = fla_lapack_slange("1", &m, &n, A, &lda, NULL);
            matrix_difference(datatype, m, n, A, lda, A_BD, lda);
            norm = fla_lapack_slange("1", &m, &n, A, &lda, NULL);
            residual = (double)(norm / (eps * norm_A * n));
            break;
        }

        case DOUBLE:
        {
            for(integer i = 0; i < nb; i++)
            {
                void *LR = NULL, *RR = NULL, *work = NULL;
                double tauQ, tauP;
                /* Collecting elementary reflector V from outputs in A */
                tauQ = ((double *)tauq)[i];
                create_vector(datatype, &LR, m);
                for(int j = 0; j < m; j++)
                {
                    if(j < i + Vt)
                        ((double *)LR)[j] = d_zero;
                    else if(j == i + Vt)
                        ((double *)LR)[j] = d_one;
                    else
                        ((double *)LR)[j] = ((double *)A_test)[j + i * lda];
                }
                /* Applying elementary reflector V to the left of the matrix */
                create_vector(datatype, &work, n);
                dlarf_("L", &m, &n, LR, &i_one, &tauQ, A, &lda, work);
                free_vector(LR);
                free_vector(work);
                /* Collecting elementary reflector U from outputs in A */
                tauP = ((double *)taup)[i];
                create_vector(datatype, &RR, n);
                for(int j = 0; j < n; j++)
                {
                    if(j < i + Ut)
                        ((double *)RR)[j] = d_zero;
                    else if(j == i + Ut)
                        ((double *)RR)[j] = d_one;
                    else
                        ((double *)RR)[j] = ((double *)A_test)[i + j * lda];
                }
                /* Applying elementary reflector U to the right of the matrix */
                create_vector(datatype, &work, m);
                dlarf_("R", &m, &n, RR, &i_one, &tauP, A, &lda, work);
                free_vector(RR);
                free_vector(work);
            }
            /* Get difference between the 2 bidiagonal forms of A */
            double norm, norm_A, eps;
            eps = fla_lapack_dlamch("P");
            norm_A = fla_lapack_dlange("1", &m, &n, A, &lda, NULL);
            matrix_difference(datatype, m, n, A, lda, A_BD, lda);
            norm = fla_lapack_dlange("1", &m, &n, A, &lda, NULL);
            residual = (double)(norm / (eps * norm_A * n));
            break;
        }

        case COMPLEX:
        {
            for(integer i = 0; i < nb; i++)
            {
                void *LR = NULL, *RR = NULL, *work = NULL;
                scomplex tauQ, tauP;
                /* Collecting elementary reflector V from outputs in A */
                tauQ = ((scomplex *)tauq)[i];
                /* Apply conjugate of TAUQ to larf to use apply V**H */
                tauQ.imag = (-1) * tauQ.imag;
                create_vector(datatype, &LR, m);
                for(int j = 0; j < m; j++)
                {
                    if(j < i + Vt)
                        ((scomplex *)LR)[j] = c_zero;
                    else if(j == i + Vt)
                        ((scomplex *)LR)[j] = c_one;
                    else
                    {
                        ((scomplex *)LR)[j].real = ((scomplex *)A_test)[j + i * lda].real;
                        ((scomplex *)LR)[j].imag = ((scomplex *)A_test)[j + i * lda].imag;
                    }
                }
                /* Applying elementary reflector V to the left of the matrix */
                create_vector(datatype, &work, n);
                clarf_("L", &m, &n, LR, &i_one, &tauQ, A, &lda, work);
                free_vector(LR);
                free_vector(work);
                /* Collecting elementary reflector U from outputs in A */
                tauP = ((scomplex *)taup)[i];
                create_vector(datatype, &RR, n);
                for(int j = 0; j < n; j++)
                {
                    if(j < i + Ut)
                        ((scomplex *)RR)[j] = c_zero;
                    else if(j == i + Ut)
                        ((scomplex *)RR)[j] = c_one;
                    else
                    {
                        ((scomplex *)RR)[j].real = ((scomplex *)A_test)[i + j * lda].real;
                        /* Output in A is of conjugate format and needs to be negated */
                        ((scomplex *)RR)[j].imag = (-1) * ((scomplex *)A_test)[i + j * lda].imag;
                    }
                }
                /* Applying elementary reflector U to the right of the matrix */
                create_vector(datatype, &work, m);
                clarf_("R", &m, &n, RR, &i_one, &tauP, A, &lda, work);
                free_vector(RR);
                free_vector(work);
            }
            /* Get difference between the 2 bidiagonal forms of A */
            float norm, norm_A, eps;
            eps = fla_lapack_slamch("P");
            norm_A = fla_lapack_clange("1", &m, &n, A, &lda, NULL);
            matrix_difference(datatype, m, n, A, lda, A_BD, lda);
            norm = fla_lapack_clange("1", &m, &n, A, &lda, NULL);
            residual = (double)(norm / (eps * norm_A * n));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            for(integer i = 0; i < nb; i++)
            {
                void *LR = NULL, *RR = NULL, *work = NULL;
                dcomplex tauQ, tauP;
                /* Collecting elementary reflector V from outputs in A */
                tauQ = ((dcomplex *)tauq)[i];
                /* Apply conjugate of TAUQ to larf to use apply V**H */
                tauQ.imag = (-1) * tauQ.imag;
                create_vector(datatype, &LR, m);
                for(int j = 0; j < m; j++)
                {
                    if(j < i + Vt)
                        ((dcomplex *)LR)[j] = z_zero;
                    else if(j == i + Vt)
                        ((dcomplex *)LR)[j] = z_one;
                    else
                    {
                        ((dcomplex *)LR)[j].real = ((dcomplex *)A_test)[j + i * lda].real;
                        ((dcomplex *)LR)[j].imag = ((dcomplex *)A_test)[j + i * lda].imag;
                    }
                }
                /* Applying elementary reflector V to the left of the matrix */
                create_vector(datatype, &work, n);
                zlarf_("L", &m, &n, LR, &i_one, &tauQ, A, &lda, work);
                free_vector(LR);
                free_vector(work);
                /* Collecting elementary reflector U from outputs in A */
                tauP = ((dcomplex *)taup)[i];
                create_vector(datatype, &RR, n);
                for(int j = 0; j < n; j++)
                {
                    if(j < i + Ut)
                        ((dcomplex *)RR)[j] = z_zero;
                    else if(j == i + Ut)
                        ((dcomplex *)RR)[j] = z_one;
                    else
                    {
                        ((dcomplex *)RR)[j].real = ((dcomplex *)A_test)[i + j * lda].real;
                        /* Output in A is of conjugate format and needs to be negated */
                        ((dcomplex *)RR)[j].imag = (-1) * ((dcomplex *)A_test)[i + j * lda].imag;
                    }
                }
                /* Applying elementary reflector U to the right of the matrix */
                create_vector(datatype, &work, m);
                zlarf_("R", &m, &n, RR, &i_one, &tauP, A, &lda, work);
                free_vector(RR);
                free_vector(work);
            }
            /* Get difference between the 2 bidiagonal forms of A */
            double norm, norm_A, eps;
            eps = fla_lapack_dlamch("P");
            norm_A = fla_lapack_zlange("1", &m, &n, A, &lda, NULL);
            matrix_difference(datatype, m, n, A, lda, A_BD, lda);
            norm = fla_lapack_zlange("1", &m, &n, A, &lda, NULL);
            residual = (double)(norm / (eps * norm_A * n));
            break;
        }
    }
    free_matrix(U);
    free_matrix(V);
    free_matrix(XU);
    free_matrix(VY);
    free_matrix(A_BD);
    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
}
void collect_LABRD_Output(integer datatype, integer m, integer n, integer nb, void *A_test,
                          integer lda, void *X, integer ldx, void *Y, integer ldy, void *U, void *V,
                          void *XU, void *VY)
{
    switch(datatype)
    {
        case FLOAT:
        {
            /* Collecting V from outputs in A */
            for(int j = 0; j < nb; j++)
            {
                for(int i = 0; i < m; i++)
                {
                    if(i < j + Vt)
                        ((float *)V)[i + j * m] = s_zero;
                    else if(i == j + Vt)
                        ((float *)V)[i + j * m] = s_one;
                    else
                        ((float *)V)[i + j * m] = ((float *)A_test)[i + j * lda];
                }
            }
            /* Collecting U**H from outputs in A */
            for(int j = 0; j < n; j++)
            {
                for(int i = 0; i < nb; i++)
                {
                    if(j < i + Ut)
                        ((float *)U)[i + j * nb] = s_zero;
                    else if(j == i + Ut)
                        ((float *)U)[i + j * nb] = s_one;
                    else
                        ((float *)U)[i + j * nb] = ((float *)A_test)[i + j * lda];
                }
            }
            /* Compute X*U**H */
            sgemm_("N", "N", &m, &n, &nb, &s_one, X, &ldx, U, &nb, &s_zero, XU, &m);
            /* Compute V*Y**H */
            sgemm_("N", "T", &m, &n, &nb, &s_one, V, &m, Y, &ldy, &s_zero, VY, &m);
            break;
        }
        case DOUBLE:
        {
            /* Collecting V from outputs in A */
            for(int j = 0; j < nb; j++)
            {
                for(int i = 0; i < m; i++)
                {
                    if(i < j + Vt)
                        ((double *)V)[i + j * m] = d_zero;
                    else if(i == j + Vt)
                        ((double *)V)[i + j * m] = d_one;
                    else
                        ((double *)V)[i + j * m] = ((double *)A_test)[i + j * lda];
                }
            }
            /* Collecting U**H from outputs in A */
            for(int j = 0; j < n; j++)
            {
                for(int i = 0; i < nb; i++)
                {
                    if(j < i + Ut)
                        ((double *)U)[i + j * nb] = d_zero;
                    else if(j == i + Ut)
                        ((double *)U)[i + j * nb] = d_one;
                    else
                        ((double *)U)[i + j * nb] = ((double *)A_test)[i + j * lda];
                }
            }
            /* Compute X*U**H */
            dgemm_("N", "N", &m, &n, &nb, &d_one, X, &ldx, U, &nb, &d_zero, XU, &m);
            /* Compute V*Y**H */
            dgemm_("N", "T", &m, &n, &nb, &d_one, V, &m, Y, &ldy, &d_zero, VY, &m);
            break;
        }
        case COMPLEX:
        {
            /* Collecting V from outputs in A */
            for(int j = 0; j < nb; j++)
            {
                for(int i = 0; i < m; i++)
                {
                    if(i < j + Vt)
                        ((scomplex *)V)[i + j * m] = c_zero;
                    else if(i == j + Vt)
                        ((scomplex *)V)[i + j * m] = c_one;
                    else
                    {
                        ((scomplex *)V)[i + j * m].real = ((scomplex *)A_test)[i + j * lda].real;
                        ((scomplex *)V)[i + j * m].imag = ((scomplex *)A_test)[i + j * lda].imag;
                    }
                }
            }
            /* Collecting U**H from outputs in A */
            for(int j = 0; j < n; j++)
            {
                for(int i = 0; i < nb; i++)
                {
                    if(j < i + Ut)
                    {
                        ((scomplex *)U)[i + j * nb] = c_zero;
                    }
                    else if(j == i + Ut)
                        ((scomplex *)U)[i + j * nb] = c_one;
                    else
                    {
                        ((scomplex *)U)[i + j * nb].real = ((scomplex *)A_test)[i + j * lda].real;
                        ((scomplex *)U)[i + j * nb].imag = ((scomplex *)A_test)[i + j * lda].imag;
                    }
                }
            }
            /* Compute X*U**H */
            cgemm_("N", "N", &m, &n, &nb, &c_one, X, &ldx, U, &nb, &c_zero, XU, &m);
            /* Compute V*Y**H */
            cgemm_("N", "C", &m, &n, &nb, &c_one, V, &m, Y, &ldy, &c_zero, VY, &m);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Collecting V from outputs in A */
            for(int j = 0; j < nb; j++)
            {
                for(int i = 0; i < m; i++)
                {
                    if(i < j + Vt)
                        ((dcomplex *)V)[i + j * m] = z_zero;
                    else if(i == j + Vt)
                        ((dcomplex *)V)[i + j * m] = z_one;
                    else
                    {
                        ((dcomplex *)V)[i + j * m].real = ((dcomplex *)A_test)[i + j * lda].real;
                        ((dcomplex *)V)[i + j * m].imag = ((dcomplex *)A_test)[i + j * lda].imag;
                    }
                }
            }
            /* Collecting U**H from outputs in A */
            for(int j = 0; j < n; j++)
            {
                for(int i = 0; i < nb; i++)
                {
                    if(j < i + Ut)
                        ((dcomplex *)U)[i + j * nb] = z_zero;
                    else if(j == i + Ut)
                        ((dcomplex *)U)[i + j * nb] = z_one;
                    else
                    {
                        ((dcomplex *)U)[i + j * nb].real = ((dcomplex *)A_test)[i + j * lda].real;
                        ((dcomplex *)U)[i + j * nb].imag = ((dcomplex *)A_test)[i + j * lda].imag;
                    }
                }
            }
            /* Compute X*U**H */
            zgemm_("N", "N", &m, &n, &nb, &z_one, X, &ldx, U, &nb, &z_zero, XU, &m);
            /* Compute V*Y**H */
            zgemm_("N", "C", &m, &n, &nb, &z_one, V, &m, Y, &ldy, &z_zero, VY, &m);
            break;
        }
    }
}
void create_bidiagonal_form(integer datatype, integer m, integer n, integer nb, void *A,
                            integer lda, void *XU, void *YV, void *d, void *e)
{
    switch(datatype)
    {
        case FLOAT:
        {
            for(integer i = 0; i < m; i++)
            {
                for(integer j = 0; j < n; j++)
                {
                    /* The first nb positions of the primary diagonal is stored in d */
                    if(i == j && i < nb)
                    {
                        ((float *)A)[i + j * lda] = ((float *)d)[i];
                    }
                    /* The first nb positions of offdiagonal is stored in e
                     * The off diagonal is present in
                     * Upper triangle of matrix if m>=n
                     * Lower triangle of matrix if m<n */
                    else if(m >= n && i + 1 == j && i < nb)
                    {
                        ((float *)A)[i + j * lda] = ((float *)e)[i];
                    }
                    else if(m < n && i - 1 == j && j < nb)
                    {
                        ((float *)A)[i + j * lda] = ((float *)e)[j];
                    }
                    /* The first nb rows and columns are reduced to 0 */
                    else if(i < nb || j < nb)
                    {
                        ((float *)A)[i + j * lda] = s_zero;
                    }
                    /* The elements of the vectors v and u together form the m-by-nb matrix
                     * V and the nb-by-n matrix U**T which are needed, with X and Y, to apply
                     * the transformation to the unreduced part of the matrix, using a block
                     * update of the form:  A := A - V*Y**T - X*U**T */
                    else
                    {
                        ((float *)A)[i + j * lda] = ((float *)A)[i + j * lda]
                                                    - ((float *)XU)[i + j * m]
                                                    - ((float *)YV)[i + j * m];
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(integer i = 0; i < m; i++)
            {
                for(integer j = 0; j < n; j++)
                {
                    /* The first nb positions of the primary diagonal is stored in d */
                    if(i == j && i < nb)
                    {
                        ((double *)A)[i + j * lda] = ((double *)d)[i];
                    }
                    /* The first nb positions of offdiagonal is stored in e
                     * The off diagonal is present in
                     * Upper triangle of matrix if m>=n
                     * Lower triangle of matrix if m<n */
                    else if(m >= n && i + 1 == j && i < nb)
                    {
                        ((double *)A)[i + j * lda] = ((double *)e)[i];
                    }
                    else if(m < n && i - 1 == j && j < nb)
                    {
                        ((double *)A)[i + j * lda] = ((double *)e)[j];
                    }
                    /* The first nb rows and columns are reduced to 0 */
                    else if(i < nb || j < nb)
                    {
                        ((double *)A)[i + j * lda] = d_zero;
                    }
                    /* The elements of the vectors v and u together form the m-by-nb matrix
                     * V and the nb-by-n matrix U**T which are needed, with X and Y, to apply
                     * the transformation to the unreduced part of the matrix, using a block
                     * update of the form:  A := A - V*Y**T - X*U**T */
                    else
                    {
                        ((double *)A)[i + j * lda] = ((double *)A)[i + j * lda]
                                                     - ((double *)XU)[i + j * m]
                                                     - ((double *)YV)[i + j * m];
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(integer i = 0; i < m; i++)
            {
                for(integer j = 0; j < n; j++)
                {
                    /* The first nb positions of the primary diagonal is stored in d */
                    if(i == j && i < nb)
                    {
                        ((scomplex *)A)[i + j * lda].real = ((float *)d)[i];
                        ((scomplex *)A)[i + j * lda].imag = s_zero;
                    }
                    /* The first nb positions of offdiagonal is stored in e
                     * The off diagonal is present in
                     * Upper triangle of matrix if m>=n
                     * Lower triangle of matrix if m<n */
                    else if(m >= n && i + 1 == j && i < nb)
                    {
                        ((scomplex *)A)[i + j * lda].real = ((float *)e)[i];
                        ((scomplex *)A)[i + j * lda].imag = s_zero;
                    }
                    else if(m < n && i - 1 == j && j < nb)
                    {
                        ((scomplex *)A)[i + j * lda].real = ((float *)e)[j];
                        ((scomplex *)A)[i + j * lda].imag = s_zero;
                    }
                    /* The first nb rows and columns are reduced to 0 */
                    else if(i < nb || j < nb)
                    {
                        ((scomplex *)A)[i + j * lda] = c_zero;
                    }
                    /* The elements of the vectors v and u together form the m-by-nb matrix
                     * V and the nb-by-n matrix U**T which are needed, with X and Y, to apply
                     * the transformation to the unreduced part of the matrix, using a block
                     * update of the form:  A := A - V*Y**T - X*U**T */
                    else
                    {
                        ((scomplex *)A)[i + j * lda].real = ((scomplex *)A)[i + j * lda].real
                                                            - ((scomplex *)XU)[i + j * m].real
                                                            - ((scomplex *)YV)[i + j * m].real;
                        ((scomplex *)A)[i + j * lda].imag = ((scomplex *)A)[i + j * lda].imag
                                                            - ((scomplex *)XU)[i + j * m].imag
                                                            - ((scomplex *)YV)[i + j * m].imag;
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(integer i = 0; i < m; i++)
            {
                for(integer j = 0; j < n; j++)
                {
                    /* The first nb positions of the primary diagonal is stored in d */
                    if(i == j && i < nb)
                    {
                        ((dcomplex *)A)[i + j * lda].real = ((double *)d)[i];
                        ((dcomplex *)A)[i + j * lda].imag = d_zero;
                    }
                    /* The first nb positions of offdiagonal is stored in e
                     * The off diagonal is present in
                     * Upper triangle of matrix if m>=n
                     * Lower triangle of matrix if m<n */
                    else if(m >= n && i + 1 == j && i < nb)
                    {
                        ((dcomplex *)A)[i + j * lda].real = ((double *)e)[i];
                        ((dcomplex *)A)[i + j * lda].imag = d_zero;
                    }
                    else if(m < n && i - 1 == j && j < nb)
                    {
                        ((dcomplex *)A)[i + j * lda].real = ((double *)e)[j];
                        ((dcomplex *)A)[i + j * lda].imag = d_zero;
                    }
                    /* The first nb rows and columns are reduced to 0 */
                    else if(i < nb || j < nb)
                    {
                        ((dcomplex *)A)[i + j * lda] = z_zero;
                    }
                    /* The elements of the vectors v and u together form the m-by-nb matrix
                     * V and the nb-by-n matrix U**T which are needed, with X and Y, to apply
                     * the transformation to the unreduced part of the matrix, using a block
                     * update of the form:  A := A - V*Y**T - X*U**T */
                    else
                    {
                        ((dcomplex *)A)[i + j * lda].real = ((dcomplex *)A)[i + j * lda].real
                                                            - ((dcomplex *)XU)[i + j * m].real
                                                            - ((dcomplex *)YV)[i + j * m].real;
                        ((dcomplex *)A)[i + j * lda].imag = ((dcomplex *)A)[i + j * lda].imag
                                                            - ((dcomplex *)XU)[i + j * m].imag
                                                            - ((dcomplex *)YV)[i + j * m].imag;
                    }
                }
            }
            break;
        }
    }
}