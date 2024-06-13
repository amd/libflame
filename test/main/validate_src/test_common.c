/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc.  All rights reserved. Portions of this
   file consist of AI-generated content.
*/

#include "test_common.h"

// Global variables
integer i_zero = 0, i_one = 1, i_n_one = -1;
float s_zero = 0, s_one = 1, s_n_one = -1;
double d_zero = 0, d_one = 1, d_n_one = -1;
scomplex c_zero = {0, 0}, c_one = {1, 0}, c_n_one = {-1, 0};
dcomplex z_zero = {0, 0}, z_one = {1, 0}, z_n_one = {-1, 0};

/* Integer absolute function */
integer fla_i_abs(integer *x)
{
    return (*x >= 0 ? (*x) : (-*x));
}
/* Allocate dynamic memory. If FLA_MEM_UNALIGNED is set, unaligned memory is allocated */
char *fla_mem_alloc(size_t size)
{
    char *buff = NULL;
#ifdef FLA_MEM_UNALIGNED
    buff = (char *)malloc(size + 1);
    if(buff == NULL)
    {
        fprintf(stderr, "malloc() returned NULL pointer\n");
        abort();
    }
    /* making aligned address to byte aligned */
    buff = buff + 1;
#else
    buff = (char *)malloc(size);
    if(buff == NULL)
    {
        fprintf(stderr, "malloc() returned NULL pointer\n");
        abort();
    }
#endif
    return buff;
}
/* create vector of given datatype*/
void create_vector(integer datatype, void **A, integer M)
{
    *A = NULL;

    switch(datatype)
    {
        case INTEGER:
        {
            *A = (integer *)fla_mem_alloc(fla_max(1, M) * sizeof(integer));
            break;
        }

        case FLOAT:
        {
            *A = (float *)fla_mem_alloc(fla_max(1, M) * sizeof(float));
            break;
        }

        case DOUBLE:
        {
            *A = (double *)fla_mem_alloc(fla_max(1, M) * sizeof(double));
            break;
        }

        case COMPLEX:
        {
            *A = (scomplex *)fla_mem_alloc(fla_max(1, M) * sizeof(scomplex));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            *A = (dcomplex *)fla_mem_alloc(fla_max(1, M) * sizeof(dcomplex));
            break;
        }
    }
}

void create_realtype_vector(integer datatype, void **A, integer M)
{
    *A = NULL;

    if(datatype == FLOAT || datatype == COMPLEX)
        *A = (float *)fla_mem_alloc(fla_max(1, M) * sizeof(float));
    else
        *A = (double *)fla_mem_alloc(fla_max(1, M) * sizeof(double));
}
/*Assign datatype*/
void assign_value(integer datatype, void *x, double data_real, double data_imag)
{
    switch(datatype)
    {
        case FLOAT:
        {
            float a = (float)data_real;
            *(float *)x = a;
            break;
        }
        case DOUBLE:
        {
            double a = data_real;
            *(double *)x = a;
            break;
        }
        case COMPLEX:
        {
            ((scomplex *)x)[0].real = (float)data_real;
            ((scomplex *)x)[0].imag = (float)data_imag;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            ((dcomplex *)x)[0].real = data_real;
            ((dcomplex *)x)[0].imag = data_imag;
            break;
        }
    }
}

/* free vector */
void free_vector(void *A)
{
    if(!A)
        return;
#ifdef FLA_MEM_UNALIGNED
    /* reset the incremented address to normal to proper freeing of memory */
    char *temp = (char *)A;
    A = (void *)(temp - 1);
#endif
    free(A);
}

/* initialize to zero */
void reset_vector(integer datatype, void *A, integer M, integer incA)
{
    integer i;

    switch(datatype)
    {
        case INTEGER:
        {
            for(i = 0; i < M; i++)
            {
                ((integer *)A)[i * incA] = 0;
            }
            break;
        }
        case FLOAT:
        {
            for(i = 0; i < M; i++)
            {
                ((float *)A)[i * incA] = 0.f;
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < M; i++)
            {
                ((double *)A)[i * incA] = 0.;
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < M; i++)
            {
                ((scomplex *)A)[i * incA].real = 0.f;
                ((scomplex *)A)[i * incA].imag = 0.f;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < M; i++)
            {
                ((dcomplex *)A)[i * incA].real = 0.;
                ((dcomplex *)A)[i * incA].imag = 0.;
            }
            break;
        }
    }
}

/* Initialize vector,
 * if range = V then initialize real type vector with random values between given range (VL, VU)
 * if range = U then initialize real type vector within specific range with uniform initialization
 * if range = R initialize vector with random values
 */
void rand_vector(integer datatype, integer M, void *A, integer inc, double VL, double VU,
                 char range)
{
    integer i;
    double step = VU - VL / M;

    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < M; i++)
            {
                if(range == 'V' || range == 'v')
                {
                    ((float *)A)[i * inc] = SRAND_IN_RANGE(VL, VU);
                }
                else if(range == 'U' || range == 'u')
                {
                    ((float *)A)[i * inc] = VL;
                    VL = VL + step;
                }
                else
                {
                    ((float *)A)[i * inc] = SRAND();
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < M; i++)
            {
                if(range == 'V' || range == 'v')
                {
                    ((double *)A)[i * inc] = DRAND_IN_RANGE(VL, VU);
                }
                else if(range == 'U' || range == 'u')
                {
                    ((double *)A)[i * inc] = VL;
                    VL = VL + step;
                }
                else
                {
                    ((double *)A)[i * inc] = DRAND();
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < M; i++)
            {
                ((scomplex *)A)[i * inc].real = SRAND();
                ((scomplex *)A)[i * inc].imag = SRAND();
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < M; i++)
            {
                ((dcomplex *)A)[i * inc].real = DRAND();
                ((dcomplex *)A)[i * inc].imag = DRAND();
            }
            break;
        }
    }
}

/* Copy a vector */
void copy_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB)
{
    switch(datatype)
    {
        case INTEGER:
        {
            integer i, iA = 0, iB = 0;

            if(LDA < 0)
            {
                iA = (-M + 1) * LDA;
            }

            if(LDB < 0)
            {
                iB = (-M + 1) * LDB;
            }

            for(i = 0; i < M; i++)
            {
                ((integer *)B)[iB] = ((integer *)A)[iA];
                iA += LDA;
                iB += LDB;
            }
            break;
        }
        case FLOAT:
        {
            scopy_(&M, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE:
        {
            dcopy_(&M, A, &LDA, B, &LDB);
            break;
        }
        case COMPLEX:
        {
            ccopy_(&M, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zcopy_(&M, A, &LDA, B, &LDB);
            break;
        }
    }
}

/* copy subvector
 * m - elements in the vector to be copied
 * A - Source matrix
 * B - Destination matrix
 * (srow, scol) - start location of the source vector in a matrix
 * (if A is a vector (srow, scol) = (0,0))
 * (drow, dcol) - start location of the destination vector in a matrix
 * (if B is a vector (drow, dcol) = (0,0)) */

void copy_subvector(integer datatype, integer m, void *A, integer lda, void *B, integer ldb,
                    integer srow, integer scol, integer drow, integer dcol)
{
    void *Ax = NULL, *Bx = NULL;
    switch(datatype)
    {
        case FLOAT:
        {
            Ax = ((float *)A + (scol * lda + srow));
            Bx = ((float *)B + (dcol * ldb + drow));
            break;
        }
        case DOUBLE:
        {
            Ax = ((double *)A + (scol * lda + srow));
            Bx = ((double *)B + (dcol * ldb + drow));
            break;
        }
        case COMPLEX:
        {
            Ax = ((scomplex *)A + (scol * lda + srow));
            Bx = ((scomplex *)B + (dcol * ldb + drow));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            Ax = ((dcomplex *)A + (scol * lda + srow));
            Bx = ((dcomplex *)B + (dcol * ldb + drow));
            break;
        }
    }
    copy_vector(datatype, m, Ax, 1, Bx, 1);
}

void copy_realtype_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB)
{
    if(datatype == FLOAT || datatype == COMPLEX)
        scopy_(&M, A, &LDA, B, &LDB);
    else
        dcopy_(&M, A, &LDA, B, &LDB);
}

/* create matrix of given datatype*/
void create_matrix(integer datatype, void **A, integer M, integer N)
{
    *A = NULL;

    switch(datatype)
    {
        case INTEGER:
        {
            *A = (integer *)fla_mem_alloc(fla_max(1, M) * fla_max(1, N) * sizeof(integer));
            break;
        }

        case FLOAT:
        {
            *A = (float *)fla_mem_alloc(fla_max(1, M) * fla_max(1, N) * sizeof(float));
            break;
        }

        case DOUBLE:
        {
            *A = (double *)fla_mem_alloc(fla_max(1, M) * fla_max(1, N) * sizeof(double));
            break;
        }

        case COMPLEX:
        {
            *A = (scomplex *)fla_mem_alloc(fla_max(1, M) * fla_max(1, N) * sizeof(scomplex));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            *A = (dcomplex *)fla_mem_alloc(fla_max(1, M) * fla_max(1, N) * sizeof(dcomplex));
            break;
        }
    }
}

void create_realtype_matrix(integer datatype, void **A, integer M, integer N)
{
    *A = NULL;

    if(datatype == FLOAT || datatype == COMPLEX)
        *A = (float *)fla_mem_alloc(fla_max(1, M) * fla_max(1, N) * sizeof(float));
    else
        *A = (double *)fla_mem_alloc(fla_max(1, M) * fla_max(1, N) * sizeof(double));
}

void *get_m_ptr(integer datatype, void *A, integer M, integer N, integer LDA)
{
    void *mat = NULL;

    switch(datatype)
    {
        case FLOAT:
        {
            mat = ((float *)A) + M + N * LDA;
            break;
        }
        case DOUBLE:
        {
            mat = ((double *)A) + M + N * LDA;
            break;
        }
        case COMPLEX:
        {
            mat = ((scomplex *)A) + M + N * LDA;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            mat = ((dcomplex *)A) + M + N * LDA;
            break;
        }
    }

    return mat;
}

/* free matrix */
void free_matrix(void *A)
{
    if(!A)
        return;
#ifdef FLA_MEM_UNALIGNED
    /* reset the incremented address to normal to proper freeing of memory */
    char *temp = (char *)A;
    A = (void *)(temp - 1);
#endif
    free(A);
}

/* Initialize matrix with random values */
void rand_matrix(integer datatype, void *A, integer M, integer N, integer LDA)
{
    integer i, j;
    if(LDA < M)
        return;
    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((float *)A)[i * LDA + j] = SRAND();
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((double *)A)[i * LDA + j] = DRAND();
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((scomplex *)A)[i * LDA + j].real = SRAND();
                    ((scomplex *)A)[i * LDA + j].imag = SRAND();
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((dcomplex *)A)[i * LDA + j].real = DRAND();
                    ((dcomplex *)A)[i * LDA + j].imag = DRAND();
                }
            }
            break;
        }
    }
}

/* Initialize symmetric matrix with random values */
void rand_sym_matrix(integer datatype, void *A, integer M, integer N, integer LDA)
{
    integer i, j;
    if(LDA < M)
        return;
    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j < M; j++)
                {
                    ((float *)A)[i * LDA + j] = SRAND();
                    ((float *)A)[j * LDA + i] = ((float *)A)[i * LDA + j];
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j < M; j++)
                {
                    ((double *)A)[i * LDA + j] = DRAND();
                    ((double *)A)[j * LDA + i] = ((double *)A)[i * LDA + j];
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j < M; j++)
                {
                    ((scomplex *)A)[i * LDA + j].real = SRAND();
                    ((scomplex *)A)[i * LDA + j].imag = SRAND();
                    ((scomplex *)A)[j * LDA + i].real = ((scomplex *)A)[i * LDA + j].real;
                    ((scomplex *)A)[j * LDA + i].imag = ((scomplex *)A)[i * LDA + j].imag;
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j < M; j++)
                {
                    ((dcomplex *)A)[i * LDA + j].real = DRAND();
                    ((dcomplex *)A)[i * LDA + j].imag = DRAND();
                    ((dcomplex *)A)[j * LDA + i].real = ((dcomplex *)A)[i * LDA + j].real;
                    ((dcomplex *)A)[j * LDA + i].imag = ((dcomplex *)A)[i * LDA + j].imag;
                }
            }
            break;
        }
    }
}

/* Copy a matrix */
void copy_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA, void *B,
                 integer LDB)
{
    if((LDA < M) || (LDB < M))
        return;

    switch(datatype)
    {
        case INTEGER:
        {
            integer i, j;

            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((integer *)B)[i * LDB + j] = ((integer *)A)[i * LDA + j];
                }
            }
            break;
        }
        case FLOAT:
        {
            fla_lapack_slacpy(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dlacpy(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_clacpy(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlacpy(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
    }
}

void copy_realtype_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA,
                          void *B, integer LDB)
{
    if(datatype == FLOAT || datatype == COMPLEX)
        fla_lapack_slacpy(uplo, &M, &N, A, &LDA, B, &LDB);
    else
        fla_lapack_dlacpy(uplo, &M, &N, A, &LDA, B, &LDB);
}

/* Initialize a matrix with zeros */
void reset_matrix(integer datatype, integer M, integer N, void *A, integer LDA)
{
    integer i, j;
    if(LDA < M)
        return;
    switch(datatype)
    {
        case INTEGER:
        {
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((integer *)A)[i * LDA + j] = 0;
                }
            }
            break;
        }

        case FLOAT:
        {
            fla_lapack_slaset("A", &M, &N, &s_zero, &s_zero, A, &LDA);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dlaset("A", &M, &N, &d_zero, &d_zero, A, &LDA);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_claset("A", &M, &N, &c_zero, &c_zero, A, &LDA);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlaset("A", &M, &N, &z_zero, &z_zero, A, &LDA);
            break;
        }
    }
}

/* Set a matrix to identity */
void set_identity_matrix(integer datatype, integer M, integer N, void *A, integer LDA)
{
    if(LDA < M)
        return;

    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slaset("A", &M, &N, &s_zero, &s_one, A, &LDA);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dlaset("A", &M, &N, &d_zero, &d_one, A, &LDA);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_claset("A", &M, &N, &c_zero, &c_one, A, &LDA);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlaset("A", &M, &N, &z_zero, &z_one, A, &LDA);
            break;
        }
    }
}

void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp)
{
    dcomplex a = *ap;
    dcomplex b = *bp;
    double temp;

    temp = b.real * b.real + b.imag * b.imag;
    if(!temp)
    {
        fprintf(stderr, "z_div_t : temp is zero. Abort\n");
        abort();
    }

    cp->real = (a.real * b.real + a.imag * b.imag) / temp;
    cp->imag = (a.imag * b.real - a.real * b.imag) / temp;
}

/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp)
{
    scomplex a = *ap;
    scomplex b = *bp;
    float temp;

    temp = b.real * b.real + b.imag * b.imag;
    if(!temp)
    {
        fprintf(stderr, "z_div_t : temp is zero. Abort\n");
        abort();
    }

    cp->real = (a.real * b.real + a.imag * b.imag) / temp;
    cp->imag = (a.imag * b.real - a.real * b.imag) / temp;
}

/* work value calculation */
integer get_work_value(integer datatype, void *work)
{
    integer value;

    if(!work)
        return 0;

    switch(datatype)
    {
        case INTEGER:
        {
            value = (*(integer *)work);
            break;
        }
        case FLOAT:
        {
            value = (integer)(*(float *)work);
            break;
        }
        case DOUBLE:
        {
            value = (integer)(*(double *)work);
            break;
        }
        case COMPLEX:
        {
            value = (integer)(((scomplex *)work)->real);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            value = (integer)(((dcomplex *)work)->real);
            break;
        }
        default:
        {
            value = 0;
            break;
        }
    }
    return value;
}

void diagmv(integer datatype, integer m, integer n, void *x, integer incx, void *a, integer a_rs,
            integer a_cs)
{
    integer inca, lda;
    integer n_iter;
    integer n_elem;
    integer j;

    if(m <= 0 || n <= 0)
        return;

    // Initialize with optimal values for column-major storage.
    inca = a_rs;
    lda = a_cs;
    n_iter = n;
    n_elem = m;

    switch(datatype)
    {
        case FLOAT:
        {
            float *a_begin;
            for(j = 0; j < n_iter; j++)
            {
                a_begin = (float *)a + j * lda;
                scalv(datatype, n_elem, x, incx, a_begin, inca);
            }
            break;
        }

        case DOUBLE:
        {
            double *a_begin;
            for(j = 0; j < n_iter; j++)
            {
                a_begin = (double *)a + j * lda;
                scalv(datatype, n_elem, x, incx, a_begin, inca);
            }
            break;
        }

        case COMPLEX:
        {
            scomplex *a_begin;
            for(j = 0; j < n_iter; j++)
            {
                a_begin = (scomplex *)a + j * lda;
                scalv(datatype, n_elem, x, incx, a_begin, inca);
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            dcomplex *a_begin;
            for(j = 0; j < n_iter; j++)
            {
                a_begin = (dcomplex *)a + j * lda;
                scalv(datatype, n_elem, x, incx, a_begin, inca);
            }
            break;
        }
    }
}

void scalv(integer datatype, integer n, void *x, integer incx, void *y, integer incy)
{
    integer i;

    switch(datatype)
    {
        case FLOAT:
        {
            float *chi, *psi;
            for(i = 0; i < n; ++i)
            {
                chi = (float *)x + i * incx;
                psi = (float *)y + i * incy;

                (*psi) = (*chi) * (*psi);
            }
            break;
        }

        case DOUBLE:
        {
            double *chi, *psi;
            for(i = 0; i < n; ++i)
            {
                chi = (double *)x + i * incx;
                psi = (double *)y + i * incy;

                (*psi) = (*chi) * (*psi);
            }
            break;
        }

        case COMPLEX:
        {
            float *chi;
            scomplex *psi;

            for(i = 0; i < n; ++i)
            {
                chi = (float *)x + i * incx;
                psi = (scomplex *)y + i * incy;

                psi->real = (*chi) * (psi)->real;
                psi->imag = (*chi) * (psi)->imag;
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double *chi;
            dcomplex *psi;
            for(i = 0; i < n; ++i)
            {
                chi = (double *)x + i * incx;
                psi = (dcomplex *)y + i * incy;

                psi->real = (*chi) * (psi)->real;
                psi->imag = (*chi) * (psi)->imag;
            }
            break;
        }
    }
}

void set_transpose(integer datatype, char *uplo, char *trans_A, char *trans_B)
{
    if(*uplo == 'L')
    {
        *trans_A = 'N';
        *trans_B = 'C';
    }
    else
    {
        *trans_A = 'C';
        *trans_B = 'N';
    }
}

void rand_spd_matrix(integer datatype, char *uplo, void **A, integer m, integer lda)
{
    void *sample = NULL;
    void *buff_A = NULL, *buff_B = NULL;
    void *a_temp = NULL;
    char trans_A, trans_B;
    if(lda < m)
        return;

    create_matrix(datatype, &sample, lda, m);
    create_matrix(datatype, &buff_A, lda, m);
    create_matrix(datatype, &buff_B, lda, m);

    reset_matrix(datatype, m, m, buff_A, lda);
    reset_matrix(datatype, m, m, buff_B, lda);

    create_matrix(datatype, &a_temp, lda, m);
    set_identity_matrix(datatype, m, m, a_temp, lda);

    /* Generate random symmetric matrix */
    rand_sym_matrix(datatype, sample, m, m, lda);

    /* Based on uplo set the transpose flag */
    set_transpose(datatype, uplo, &trans_A, &trans_B);

    copy_matrix(datatype, uplo, m, m, sample, lda, buff_A, lda);
    copy_matrix(datatype, uplo, m, m, sample, lda, buff_B, lda);

    switch(datatype)
    {
        case FLOAT:
        {
            float beta = m;
            sgemm_(&trans_A, &trans_B, &m, &m, &m, &s_one, buff_A, &lda, buff_B, &lda, &beta,
                   a_temp, &lda);
            break;
        }
        case DOUBLE:
        {
            double beta = m;
            dgemm_(&trans_A, &trans_B, &m, &m, &m, &d_one, buff_A, &lda, buff_B, &lda, &beta,
                   a_temp, &lda);
            break;
        }
        case COMPLEX:
        {
            scomplex beta = {m, 0};
            cgemm_(&trans_A, &trans_B, &m, &m, &m, &c_one, buff_A, &lda, buff_B, &lda, &beta,
                   a_temp, &lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex beta = {m, 0};
            zgemm_(&trans_A, &trans_B, &m, &m, &m, &z_one, buff_A, &lda, buff_B, &lda, &beta,
                   a_temp, &lda);
            break;
        }
    }
    copy_matrix(datatype, "full", m, m, a_temp, lda, *A, lda);

    /* free buffers */
    free_matrix(sample);
    free_matrix(buff_A);
    free_matrix(buff_B);
    free_matrix(a_temp);
}

/* Create diagonal matrix by copying elements from a realtype vector to matrix */
void diagonalize_realtype_vector(integer datatype, void *s, void *sigma, integer m, integer n,
                                 integer LDA)
{
    integer incr, i, j, min_m_n;

    incr = m + 1;
    min_m_n = fla_min(m, n);

    reset_matrix(datatype, m, n, sigma, m);

    switch(datatype)
    {
        case FLOAT:
        {
            scopy_(&min_m_n, s, &i_one, sigma, &incr);
            break;
        }
        case DOUBLE:
        {
            dcopy_(&min_m_n, s, &i_one, sigma, &incr);
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < n; i++)
            {
                for(j = i; j < m; j++)
                {
                    if(i == j)
                        ((scomplex *)sigma)[i * LDA + j].real = ((float *)s)[i];
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < n; i++)
            {
                for(j = i; j < m; j++)
                {
                    if(i == j)
                        ((dcomplex *)sigma)[i * LDA + j].real = ((double *)s)[i];
                }
            }
            break;
        }
    }
}

/* Generate random Hermitian matrix */
void rand_hermitian_matrix(integer datatype, integer n, void **A, integer lda)
{
    void *B = NULL;
    if(lda < n)
        return;

    create_matrix(datatype, &B, n, n);
    reset_matrix(datatype, n, n, B, n);
    rand_matrix(datatype, B, n, n, n);

    switch(datatype)
    {
        case COMPLEX:
        {
            cgemm_("N", "C", &n, &n, &n, &c_one, B, &n, B, &n, &c_zero, *A, &lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zgemm_("N", "C", &n, &n, &n, &z_one, B, &n, B, &n, &z_zero, *A, &lda);
            break;
        }
    }
    free_matrix(B);
}

/* block diagonal matrix is required for computing eigen decomposition of non symmetric matrix.
   W is a block diagonal matrix, with a 1x1 block for each
   real eigenvalue and a 2x2 block for each complex conjugate
   pair.then the 2 x 2 block corresponding to the pair will be:

              (  wr  wi  )
              ( -wi  wr  )
*/
void create_block_diagonal_matrix(integer datatype, void *wr, void *wi, void *lambda, integer m,
                                  integer n, integer lda)
{
    integer i, j;

    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < m; i++)
            {
                j = i;
                if(((float *)wi)[i] != 0.f)
                {
                    ((float *)lambda)[i * lda + j] = ((float *)wr)[i];
                    ((float *)lambda)[i * lda + (j + 1)] = -((float *)wi)[i];
                    i++;
                    j++;
                    ((float *)lambda)[i * lda + j] = ((float *)wr)[i];
                    ((float *)lambda)[i * lda + (j - 1)] = -((float *)wi)[i];
                }
                else
                {
                    ((float *)lambda)[i * lda + j] = ((float *)wr)[i];
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < m; i++)
            {
                j = i;
                if(((double *)wi)[i] != 0.)
                {
                    ((double *)lambda)[i * lda + j] = ((double *)wr)[i];
                    ((double *)lambda)[i * lda + (j + 1)] = -((double *)wi)[i];
                    i++;
                    j++;
                    ((double *)lambda)[i * lda + j] = ((double *)wr)[i];
                    ((double *)lambda)[i * lda + (j - 1)] = -((double *)wi)[i];
                }
                else
                {
                    ((double *)lambda)[i * lda + j] = ((double *)wr)[i];
                }
            }
            break;
        }
    }
}

/* If trn == 'N' => checks whether A * A**T == I
   If trn == 'T' => checks whether A**T * A == I*/

double check_orthogonal_matrix(char trn, integer datatype, void *A, integer m, integer n, integer k,
                               integer lda)
{
    void *a_temp = NULL, *work = NULL;
    double resid = 0.;
    /* Create Identity matrix to validate orthogonal property of matrix A*/
    create_matrix(datatype, &a_temp, k, k);

    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm;
            eps = fla_lapack_slamch("P");

            fla_lapack_slaset("full", &k, &k, &s_zero, &s_one, a_temp, &k);
            if(trn == 'N')
            {
                sgemm_("N", "T", &m, &m, &n, &s_one, A, &lda, A, &lda, &s_n_one, a_temp, &k);
            }
            else if(trn == 'T')
            {
                sgemm_("T", "N", &m, &m, &n, &s_one, A, &lda, A, &lda, &s_n_one, a_temp, &k);
            }
            norm = fla_lapack_slange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * fla_max(m, n)));
            break;
        }
        case DOUBLE:
        {
            double eps, norm;
            eps = fla_lapack_dlamch("P");

            fla_lapack_dlaset("full", &k, &k, &d_zero, &d_one, a_temp, &k);
            if(trn == 'N')
            {
                dgemm_("N", "T", &m, &m, &n, &d_one, A, &lda, A, &lda, &d_n_one, a_temp, &k);
            }
            else if(trn == 'T')
            {
                dgemm_("T", "N", &m, &m, &n, &d_one, A, &lda, A, &lda, &d_n_one, a_temp, &k);
            }
            norm = fla_lapack_dlange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * fla_max(m, n)));
            break;
        }
        case COMPLEX:
        {
            float eps, norm;
            eps = fla_lapack_slamch("P");
            fla_lapack_claset("full", &k, &k, &c_zero, &c_one, a_temp, &k);
            if(trn == 'N')
            {
                cgemm_("N", "C", &m, &m, &n, &c_one, A, &lda, A, &lda, &c_n_one, a_temp, &k);
            }
            else if(trn == 'C')
            {
                cgemm_("C", "N", &m, &m, &n, &c_one, A, &lda, A, &lda, &c_n_one, a_temp, &k);
            }
            norm = fla_lapack_clange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * fla_max(m, n)));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm;
            eps = fla_lapack_dlamch("P");
            fla_lapack_zlaset("full", &k, &k, &z_zero, &z_one, a_temp, &k);
            if(trn == 'N')
            {
                zgemm_("N", "C", &m, &m, &n, &z_one, A, &lda, A, &lda, &z_n_one, a_temp, &k);
            }
            else if(trn == 'C')
            {
                zgemm_("C", "N", &m, &m, &n, &z_one, A, &lda, A, &lda, &z_n_one, a_temp, &k);
            }
            norm = fla_lapack_zlange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * fla_max(m, n)));
            break;
        }
    }
    free_matrix(a_temp);
    return resid;
}

/* Checks whether A**T * A == I */
double check_orthogonality(integer datatype, void *A, integer m, integer n, integer lda)
{
    void *a_temp = NULL, *work = NULL;
    double resid = 0.;
    integer k;

    /* Create Identity matrix to validate orthogonal property of matrix A*/
    if(m <= n)
    {
        create_matrix(datatype, &a_temp, m, m);
        k = m;
    }
    else
    {
        create_matrix(datatype, &a_temp, n, n);
        k = n;
    }
    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm;
            eps = fla_lapack_slamch("P");

            fla_lapack_slaset("full", &k, &k, &s_zero, &s_one, a_temp, &k);
            sgemm_("T", "N", &k, &k, &m, &s_one, A, &lda, A, &lda, &s_n_one, a_temp, &k);
            norm = fla_lapack_slange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * (float)k));
            break;
        }
        case DOUBLE:
        {
            double eps, norm;
            eps = fla_lapack_dlamch("P");

            fla_lapack_dlaset("full", &k, &k, &d_zero, &d_one, a_temp, &k);
            dgemm_("T", "N", &k, &k, &m, &d_one, A, &lda, A, &lda, &d_n_one, a_temp, &k);
            norm = fla_lapack_dlange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * (float)k));
            break;
        }
        case COMPLEX:
        {
            float eps, norm;
            eps = fla_lapack_slamch("P");

            fla_lapack_claset("full", &k, &k, &c_zero, &c_one, a_temp, &k);
            cgemm_("C", "N", &k, &k, &m, &c_one, A, &lda, A, &lda, &c_n_one, a_temp, &k);
            norm = fla_lapack_clange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * (float)k));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm;
            eps = fla_lapack_dlamch("P");

            fla_lapack_zlaset("full", &k, &k, &z_zero, &z_one, a_temp, &k);
            zgemm_("C", "N", &k, &k, &m, &z_one, A, &lda, A, &lda, &z_n_one, a_temp, &k);
            norm = fla_lapack_zlange("1", &k, &k, a_temp, &k, work);
            resid = (double)(norm / (eps * (float)k));
            break;
        }
    }
    free_matrix(a_temp);
    return resid;
}

/* copy submatrix from a matrix
 * (m, n) - dimensions of the sub-matrix to be copied
 * A - original matirx
 * B - destination matrix
 * (srow, scol) - start location of the original matrix from where the value has to be copied
 * (drow, dcol) - start location of the destination matrix to where the values has to be copied*/

void copy_submatrix(integer datatype, integer m, integer n, void *A, integer lda, void *B,
                    integer ldb, integer srow, integer scol, integer drow, integer dcol)
{
    void *sub_A = NULL, *sub_B = NULL;
    switch(datatype)
    {
        case FLOAT:
        {
            sub_A = ((float *)A + (scol * lda + srow));
            sub_B = ((float *)B + (dcol * ldb + drow));
            break;
        }
        case DOUBLE:
        {
            sub_A = ((double *)A + (scol * lda + srow));
            sub_B = ((double *)B + (dcol * ldb + drow));
            break;
        }
        case COMPLEX:
        {
            sub_A = ((scomplex *)A + (scol * lda + srow));
            sub_B = ((scomplex *)B + (dcol * ldb + drow));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            sub_A = ((dcomplex *)A + (scol * lda + srow));
            sub_B = ((dcomplex *)B + (dcol * ldb + drow));
            break;
        }
    }
    copy_matrix(datatype, "FULL", m, n, sub_A, lda, sub_B, ldb);
}

void scgemv(char TRANS, integer real_alpha, integer m, integer n, scomplex *alpha, float *a,
            integer lda, scomplex *v, integer incv, float beta, scomplex *c, integer inc)
{
    integer i, j;
    float real, imag;
    float rl, ig;
    float alphar;
    void *A = NULL;

    create_matrix(FLOAT, &A, lda, n);

    if(TRANS == 'T')
    {
        /* Transpose of a matrix A */
        for(i = 0; i < n; i++)
        {
            for(j = 0; j < n; j++)
            {
                ((float *)A)[i * lda + j] = a[i + j * lda];
            }
        }
    }
    else
    {
        copy_matrix(FLOAT, "full", n, n, a, lda, A, lda);
    }

    if(real_alpha)
    {
        alphar = alpha->real;
        for(i = 0; i < m; i++)
        {
            real = 0;
            imag = 0;
            for(j = 0; j < n; j++)
            {
                real = real + ((float *)A)[i + j * lda] * v[j * incv].real;
                imag = imag + ((float *)A)[i + j * lda] * v[j * incv].imag;
            }
            c[i * inc].real = alphar * real + beta * c[i * inc].real;
            c[i * inc].imag = alphar * imag + beta * c[i * inc].imag;
        }
    }
    else
    {
        for(i = 0; i < m; i++)
        {
            real = 0;
            imag = 0;
            for(j = 0; j < n; j++)
            {
                real = real + ((float *)A)[i + j * lda] * v[j * incv].real;
                imag = imag + ((float *)A)[i + j * lda] * v[j * incv].imag;
            }

            rl = alpha->real * real - alpha->imag * imag;
            ig = alpha->real * imag + alpha->imag * real;

            c[i * inc].real = rl + beta * c[i * inc].real;
            c[i * inc].imag = ig + beta * c[i * inc].imag;
        }
    }

    free_matrix(A);
}

/* Get datatype from string. */
integer get_datatype(char stype)
{
    integer datatype;
    if(stype == 's' || stype == 'S')
        datatype = FLOAT;
    else if(stype == 'd' || stype == 'D')
        datatype = DOUBLE;
    else if(stype == 'c' || stype == 'C')
        datatype = COMPLEX;
    else if(stype == 'z' || stype == 'Z')
        datatype = DOUBLE_COMPLEX;
    else
        datatype = INVALID_TYPE;

    return datatype;
}

/* Get realtype of given datatype. */
integer get_realtype(integer datatype)
{
    if(datatype == FLOAT || datatype == COMPLEX)
    {
        return FLOAT;
    }
    else if(datatype == DOUBLE || datatype == DOUBLE_COMPLEX)
    {
        return DOUBLE;
    }
    else
    {
        fprintf(stderr, "Invalid datatype is passed.\n");
        return -1;
    }
}

/* Get value from pointer with realtype of given datatype. */
double get_realtype_value(integer datatype, void *value)
{
    if(datatype == FLOAT || datatype == COMPLEX)
    {
        return *(float *)value;
    }
    else if(datatype == DOUBLE || datatype == DOUBLE_COMPLEX)
    {
        return *(double *)value;
    }
    else
    {
        fprintf(stderr, "Invalid datatype is passed, returning 0.\n");
        return 0;
    }
}

/* Initialize symmetric tridiagonal matrix with random values.
   Initializes random values only for diagonal and off diagonal elements.*/
void rand_sym_tridiag_matrix(integer datatype, void *A, integer M, integer N, integer LDA)
{
    integer i, j;
    if(LDA < M)
        return;
    reset_matrix(datatype, M, N, A, LDA);

    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1; j++)
                {
                    if(j < N)
                    {
                        ((float *)A)[i * LDA + j] = SRAND();
                        ((float *)A)[j * LDA + i] = ((float *)A)[i * LDA + j];
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1; j++)
                {
                    if(j < N)
                    {
                        ((double *)A)[i * LDA + j] = DRAND();
                        ((double *)A)[j * LDA + i] = ((double *)A)[i * LDA + j];
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1; j++)
                {
                    if(j < N)
                    {
                        ((scomplex *)A)[i * LDA + j].real = SRAND();
                        ((scomplex *)A)[i * LDA + j].imag = SRAND();
                        ((scomplex *)A)[j * LDA + i].real = ((scomplex *)A)[i * LDA + j].real;
                        ((scomplex *)A)[j * LDA + i].imag = ((scomplex *)A)[i * LDA + j].imag;
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1; j++)
                {
                    if(j < N)
                    {
                        ((dcomplex *)A)[i * LDA + j].real = DRAND();
                        ((dcomplex *)A)[i * LDA + j].imag = DRAND();
                        ((dcomplex *)A)[j * LDA + i].real = ((dcomplex *)A)[i * LDA + j].real;
                        ((dcomplex *)A)[j * LDA + i].imag = ((dcomplex *)A)[i * LDA + j].imag;
                    }
                }
            }
            break;
        }
    }
}

/* Get diagonal elements of matrix A into Diag vector. */
void get_diagonal(integer datatype, void *A, integer m, integer n, integer lda, void *Diag)
{
    integer i, j;
    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0, j = 0; i < m; i++, j++)
            {
                ((float *)Diag)[i] = ((float *)A)[i * lda + j];
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0, j = 0; i < m; i++, j++)
            {
                ((double *)Diag)[i] = ((double *)A)[i * lda + j];
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0, j = 0; i < m; i++, j++)
            {
                ((scomplex *)Diag)[i].real = ((scomplex *)A)[i * lda + j].real;
                ((scomplex *)Diag)[i].imag = ((scomplex *)A)[i * lda + j].imag;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0, j = 0; i < m; i++, j++)
            {
                ((dcomplex *)Diag)[i].real = ((dcomplex *)A)[i * lda + j].real;
                ((dcomplex *)Diag)[i].imag = ((dcomplex *)A)[i * lda + j].imag;
            }
            break;
        }
    }
}

/* Get subdiagonal elements of matrix A into Subdiag vector.*/
void get_subdiagonal(integer datatype, void *A, integer m, integer n, integer lda, void *Subdiag)
{
    integer i, j;
    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 1, j = 0; i < m; i++, j++)
            {
                ((float *)Subdiag)[j] = ((float *)A)[i * lda + j];
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 1, j = 0; i < m; i++, j++)
            {
                ((double *)Subdiag)[j] = ((double *)A)[i * lda + j];
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 1, j = 0; i < m; i++, j++)
            {
                ((scomplex *)Subdiag)[j].real = ((scomplex *)A)[i * lda + j].real;
                ((scomplex *)Subdiag)[j].imag = ((scomplex *)A)[i * lda + j].imag;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 1, j = 0; i < m; i++, j++)
            {
                ((dcomplex *)Subdiag)[j].real = ((dcomplex *)A)[i * lda + j].real;
                ((dcomplex *)Subdiag)[j].imag = ((dcomplex *)A)[i * lda + j].imag;
            }
            break;
        }
    }
}

void copy_sym_tridiag_matrix(integer datatype, void *D, void *E, integer M, integer N, void *B,
                             integer LDA)
{
    integer i, j;
    if(LDA < M)
        return;
    reset_matrix(datatype, M, N, B, LDA);

    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1 && j < N; j++)
                {
                    if(j == i)
                    {
                        ((float *)B)[i * LDA + j] = ((float *)D)[i];
                    }
                    else
                    {
                        ((float *)B)[i * LDA + j] = ((float *)E)[i];
                        ((float *)B)[j * LDA + i] = ((float *)E)[i];
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1 && j < N; j++)
                {
                    if(j == i)
                    {
                        ((double *)B)[i * LDA + j] = ((double *)D)[i];
                    }
                    else
                    {
                        ((double *)B)[i * LDA + j] = ((double *)E)[i];
                        ((double *)B)[j * LDA + i] = ((double *)E)[i];
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1 && j < N; j++)
                {
                    if(j == i)
                    {
                        ((scomplex *)B)[i * LDA + j].real = ((float *)D)[i];
                    }
                    else
                    {
                        ((scomplex *)B)[i * LDA + j].real = ((float *)E)[i];
                        ((scomplex *)B)[j * LDA + i].real = ((float *)E)[i];
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < N; i++)
            {
                for(j = i; j <= i + 1 && j < N; j++)
                {
                    if(j == i)
                    {
                        ((dcomplex *)B)[i * LDA + j].real = ((double *)D)[i];
                    }
                    else
                    {
                        ((dcomplex *)B)[i * LDA + j].real = ((double *)E)[i];
                        ((dcomplex *)B)[j * LDA + i].real = ((double *)E)[i];
                    }
                }
            }
            break;
        }
    }
}
/* Get the maximum value from the array */
void get_max_from_array(integer datatype, void *arr, void *max_val, integer n)
{
    integer i;

    switch(datatype)
    {
        case INTEGER:
        {
            integer *ptr = arr;
            integer maxlocal = INT_MIN;
            for(i = 0; i < n; i++)
            {
                int val = fla_i_abs(&(ptr[i]));
                if(val > maxlocal)
                {
                    maxlocal = val;
                }
            }
            *(integer *)max_val = maxlocal;
            break;
        }

        case FLOAT:
        {
            float *ptr = arr;
            float maxlocal = FLT_MIN;
            for(i = 0; i < n; i++)
            {
                float val = FLA_FABS(ptr[i]);
                if(val > maxlocal)
                {
                    maxlocal = val;
                }
            }

            *(float *)max_val = maxlocal;
            break;
        }

        case DOUBLE:
        {
            double *ptr = arr;
            double maxlocal = DBL_MIN;
            for(i = 0; i < n; i++)
            {
                double val = FLA_FABS(ptr[i]);
                if(val > maxlocal)
                {
                    maxlocal = val;
                }
            }

            *(double *)max_val = maxlocal;
            break;
        }

        /* Implementation of complex needs to be relook*/
        case COMPLEX:
        {
            scomplex *ptr = arr;
            scomplex maxlocal;
            maxlocal.real = FLT_MIN;
            maxlocal.imag = FLT_MIN;

            for(i = 0; i < n; i++)
            {
                float real, imag;
                real = FLA_FABS(ptr[i].real);
                imag = FLA_FABS(ptr[i].imag);
                /* Compare real part */
                if(real != s_zero && real > maxlocal.real)
                {
                    maxlocal.real = real;
                }
                /* Compare imaginary part */
                if(imag != s_zero && imag > maxlocal.imag)
                {
                    maxlocal.imag = imag;
                }
            }

            *(float *)max_val = fla_max(maxlocal.real, maxlocal.imag);
            break;
        }

        /* Implementation of complex needs to be relook*/
        case DOUBLE_COMPLEX:
        {
            dcomplex *ptr = arr;
            dcomplex maxlocal;
            maxlocal.real = DBL_MIN;
            maxlocal.imag = DBL_MIN;

            for(i = 0; i < n; i++)
            {
                double real, imag;
                real = FLA_FABS(ptr[i].real);
                imag = FLA_FABS(ptr[i].imag);
                /* Compare real part */
                if(real != d_zero && real > maxlocal.real)
                {
                    maxlocal.real = real;
                }
                /* Compare imaginary part */
                if(imag != d_zero && imag > maxlocal.imag)
                {
                    maxlocal.imag = imag;
                }
            }

            *(double *)max_val = fla_max(maxlocal.real, maxlocal.imag);
            break;
        }
    }
}

/* Get the minimum value from the array */
void get_min_from_array(integer datatype, void *arr, void *min_val, integer n)
{
    integer i;

    switch(datatype)
    {
        case INTEGER:
        {
            integer *ptr = arr;
            integer minlocal = INT_MAX;

            for(i = 0; i < n; i++)
            {
                int val = fla_i_abs(&(ptr[i]));
                if(val < minlocal)
                {
                    minlocal = val;
                }
            }

            *(integer *)min_val = minlocal;
            break;
        }

        case FLOAT:
        {
            float *ptr = arr;
            float minlocal = FLT_MAX;

            for(i = 0; i < n; i++)
            {
                float val = FLA_FABS(ptr[i]);

                if(val != s_zero && val < minlocal)
                {
                    minlocal = val;
                }
            }

            *(float *)min_val = minlocal;
            break;
        }

        case DOUBLE:
        {
            double *ptr = arr;
            double minlocal = DBL_MAX;
            for(i = 0; i < n; i++)
            {
                double val = FLA_FABS(ptr[i]);
                if(val != d_zero && val < minlocal)
                {
                    minlocal = val;
                }
            }

            *(double *)min_val = minlocal;
            break;
        }

        /* Implementation of complex needs to be relook*/
        case COMPLEX:
        {
            scomplex *ptr = arr;
            scomplex minlocal;
            minlocal.real = FLT_MAX;
            minlocal.imag = FLT_MAX;

            for(i = 0; i < n; i++)
            {
                float real, imag;
                real = FLA_FABS(ptr[i].real);
                imag = FLA_FABS(ptr[i].imag);
                /* Compare real part */
                if(real != s_zero && real < minlocal.real)
                {
                    minlocal.real = real;
                }
                /* Compare imaginary part */
                if(imag != s_zero && imag < minlocal.imag)
                {
                    minlocal.imag = imag;
                }
            }

            *(float *)min_val = fla_min(minlocal.real, minlocal.imag);
            break;
        }

        /* Implementation of complex needs to be relook*/
        case DOUBLE_COMPLEX:
        {
            dcomplex *ptr = arr;
            dcomplex minlocal;
            minlocal.real = DBL_MAX;
            minlocal.imag = DBL_MAX;

            for(i = 0; i < n; i++)
            {
                double real, imag;
                real = FLA_FABS(ptr[i].real);
                imag = FLA_FABS(ptr[i].imag);
                /* Compare real part */
                if(real != d_zero && real < minlocal.real)
                {
                    minlocal.real = real;
                }
                /* Compare imaginary part */
                if(imag != d_zero && imag < minlocal.imag)
                {
                    minlocal.imag = imag;
                }
            }

            *(double *)min_val = fla_min(minlocal.real, minlocal.imag);
            break;
        }
    }
}
/* Reading matrix input data from a file */
void init_matrix_from_file(integer datatype, void *A, integer m, integer n, integer lda, FILE *fptr)
{
    int i, j;
    if(lda < m)
        return;
    switch(datatype)
    {
        case FLOAT:
        {
            float num;

            for(i = 0; i < m; i++)
            {
                for(j = 0; j < n; j++)
                {
                    fscanf(fptr, "%f", &num);
                    ((float *)A)[i * lda + j] = num;
                }
            }
            break;
        }

        case DOUBLE:
        {
            double num;

            for(i = 0; i < m; i++)
            {
                for(j = 0; j < n; j++)
                {
                    fscanf(fptr, "%lf", &num);
                    ((double *)A)[i * lda + j] = num;
                }
            }
            break;
        }

        case COMPLEX:
        {
            float num;

            for(i = 0; i < m; i++)
            {
                for(j = 0; j < n; j++)
                {
                    fscanf(fptr, "%f", &num);
                    ((scomplex *)A)[i * lda + j].real = num;
                    fscanf(fptr, "%f", &num);
                    ((scomplex *)A)[i * lda + j].imag = num;
                }
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double num;

            for(i = 0; i < m; i++)
            {
                for(j = 0; j < n; j++)
                {
                    fscanf(fptr, "%lf", &num);
                    ((dcomplex *)A)[i * lda + j].real = num;
                    fscanf(fptr, "%lf", &num);
                    ((dcomplex *)A)[i * lda + j].imag = num;
                }
            }
            break;
        }
    }
}
/* Reading vector input data from a file */
void init_vector_from_file(integer datatype, void *A, integer m, integer inc, FILE *fptr)
{
    int i;

    switch(datatype)
    {
        case FLOAT:
        {
            float num;

            for(i = 0; i < m; i++)
            {
                fscanf(fptr, "%f", &num);
                ((float *)A)[i * inc] = num;
            }
            break;
        }

        case DOUBLE:
        {
            double num;

            for(i = 0; i < m; i++)
            {
                fscanf(fptr, "%lf", &num);
                ((double *)A)[i * inc] = num;
            }
            break;
        }

        case COMPLEX:
        {
            float num;

            for(i = 0; i < m; i++)
            {
                fscanf(fptr, "%f", &num);
                ((scomplex *)A)[i * inc].real = num;
                fscanf(fptr, "%f", &num);
                ((scomplex *)A)[i * inc].imag = num;
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double num;

            for(i = 0; i < m; i++)
            {
                fscanf(fptr, "%lf", &num);
                ((dcomplex *)A)[i * inc].real = num;
                fscanf(fptr, "%lf", &num);
                ((dcomplex *)A)[i * inc].imag = num;
            }
            break;
        }
    }
}

/* Convert matrix according to ILO and IHI values */
void get_generic_triangular_matrix(integer datatype, integer N, void *A, integer LDA, integer ilo,
                                   integer ihi)
{
    if(LDA < N)
        return;
    /* Intialize matrix with random values */
    rand_matrix(datatype, A, N, N, LDA);
    integer i;

    switch(datatype)
    {
        case FLOAT:
        {
            /* Making elements below diagonal for columns 0 to ilo-1 to Zero*/
            for(i = 0; i < ilo - 1; i++)
            {
                float *p = &((float *)A)[(i + 1) + i * LDA];
                reset_vector(datatype, (void *)p, N - i - 1, 1);
            }
            /* Making elements below diagonal for rows ihi+1 to N to Zero*/
            for(i = ihi; i < N; i++)
            {
                float *p = &((float *)A)[i];
                reset_vector(datatype, (void *)p, i, LDA);
            }
            break;
        }
        case DOUBLE:
        {
            /* Making elements below diagonal for columns 0 to ilo-1 to Zero*/
            for(i = 0; i < ilo - 1; i++)
            {
                double *p = &((double *)A)[(i + 1) + i * LDA];
                reset_vector(datatype, (void *)p, N - i - 1, 1);
            }
            /* Making elements below diagonal for rows ihi+1 to N to Zero*/
            for(i = ihi; i < N; i++)
            {
                double *p = &((double *)A)[i];
                reset_vector(datatype, (void *)p, i, LDA);
            }
            break;
        }
        case COMPLEX:
        {
            /* Making elements below diagonal for columns 0 to ilo-1 to Zero*/
            for(i = 0; i < ilo - 1; i++)
            {
                scomplex *p = &((scomplex *)A)[(i + 1) + i * LDA];
                reset_vector(datatype, (void *)p, N - i - 1, 1);
            }
            /* Making elements below diagonal for rows ihi+1 to N to Zero*/
            for(i = ihi; i < N; i++)
            {
                scomplex *p = &((scomplex *)A)[i];
                reset_vector(datatype, (void *)p, i, LDA);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Making elements below diagonal for columns 0 to ilo-1 to Zero*/
            for(i = 0; i < ilo - 1; i++)
            {
                dcomplex *p = &((dcomplex *)A)[(i + 1) + i * LDA];
                reset_vector(datatype, (void *)p, N - i - 1, 1);
            }
            /* Making elements below diagonal for rows ihi+1 to N to Zero*/
            for(i = ihi; i < N; i++)
            {
                dcomplex *p = &((dcomplex *)A)[i];
                reset_vector(datatype, (void *)p, i, LDA);
            }
            break;
        }
    }
}

/* Generate Hessenberg matrix */
void get_hessenberg_matrix(integer datatype, integer n, void *A, integer lda, void *Z, integer ldz,
                           integer *ilo, integer *ihi, void *scale, integer *info)
{
    static integer g_lwork;
    void *A_save = NULL;
    void *tau = NULL, *work = NULL;
    integer lwork;
    if((lda < n) || (ldz < n))
        return;
    create_matrix(datatype, &A_save, lda, n);
    create_vector(datatype, &tau, n - 1);

    /* Convert matrix according to ILO and IHI values */
    get_generic_triangular_matrix(datatype, n, A, lda, *ilo, *ihi);

    switch(datatype)
    {
        case FLOAT:
        {
            /* Make a workspace query the first time through. This will provide us with
                    and ideal workspace size based on an internal block size.*/
            if(g_lwork <= 0)
            {
                create_vector(datatype, &work, 1);
                lwork = -1;
                fla_lapack_sgehrd(&n, ilo, ihi, NULL, &lda, NULL, work, &lwork, info);
                if(*info == 0)
                {
                    lwork = get_work_value(datatype, work);
                    free_vector(work);
                }
            }
            else
            {
                lwork = g_lwork;
            }
            create_vector(datatype, &work, lwork);

            /* Call to SGEHRD API to generate hessenberg matrix*/
            fla_lapack_sgehrd(&n, ilo, ihi, A, &lda, tau, work, &lwork, info);
            reset_vector(datatype, work, lwork, 1);
            copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

            /* Call to SORGHR API to generate orthogonal matrix*/
            fla_lapack_sorghr(&n, ilo, ihi, A_save, &lda, tau, work, &lwork, info);
            copy_matrix(datatype, "full", n, n, A_save, lda, Z, ldz);

            /* Convert matrix from SGEHRD to upper hessenberg matrix */
            convert_upper_hessenberg(datatype, n, A, lda);

            free_vector(work);
            break;
        }
        case DOUBLE:
        {
            /* Make a workspace query the first time through. This will provide us with
                    and ideal workspace size based on an internal block size.*/
            if(g_lwork <= 0)
            {
                create_vector(datatype, &work, 1);
                lwork = -1;
                fla_lapack_dgehrd(&n, ilo, ihi, NULL, &lda, NULL, work, &lwork, info);
                if(*info == 0)
                {
                    lwork = get_work_value(datatype, work);
                    free_vector(work);
                }
            }
            else
            {
                lwork = g_lwork;
            }

            create_vector(datatype, &work, lwork);

            /* Call to DGEHRD API to generate hessenberg matrix*/
            fla_lapack_dgehrd(&n, ilo, ihi, A, &lda, tau, work, &lwork, info);
            reset_vector(datatype, work, lwork, 1);
            copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

            /* Call to DORGHR API to generate orthogonal matrix*/
            fla_lapack_dorghr(&n, ilo, ihi, A_save, &lda, tau, work, &lwork, info);
            copy_matrix(datatype, "full", n, n, A_save, lda, Z, ldz);

            /* Convert matrix from DGEHRD to upper hessenberg matrix */
            convert_upper_hessenberg(datatype, n, A, lda);

            free_vector(work);
            break;
        }
        case COMPLEX:
        {
            /* Make a workspace query the first time through. This will provide us with
                    and ideal workspace size based on an internal block size.*/
            if(g_lwork <= 0)
            {
                create_vector(datatype, &work, 1);
                lwork = -1;
                fla_lapack_cgehrd(&n, ilo, ihi, NULL, &lda, NULL, work, &lwork, info);
                if(*info == 0)
                {
                    lwork = get_work_value(datatype, work);
                    free_vector(work);
                }
            }
            else
            {
                lwork = g_lwork;
            }
            create_vector(datatype, &work, lwork);

            /* Call to CGEHRD API to generate hessenberg matrix*/
            fla_lapack_cgehrd(&n, ilo, ihi, A, &lda, tau, work, &lwork, info);
            reset_vector(datatype, work, lwork, 1);
            copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

            /* Call to CUNGHR API to generate orthogonal matrix*/
            fla_lapack_cunghr(&n, ilo, ihi, A_save, &lda, tau, work, &lwork, info);
            copy_matrix(datatype, "full", n, n, A_save, lda, Z, ldz);

            /* Convert matrix from CGEHRD to upper hessenberg matrix */
            convert_upper_hessenberg(datatype, n, A, lda);

            free_vector(work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Make a workspace query the first time through. This will provide us with
                    and ideal workspace size based on an internal block size.*/
            if(g_lwork <= 0)
            {
                create_vector(datatype, &work, 1);
                lwork = -1;
                fla_lapack_zgehrd(&n, ilo, ihi, NULL, &lda, NULL, work, &lwork, info);
                if(*info == 0)
                {
                    lwork = get_work_value(datatype, work);
                    free_vector(work);
                }
            }
            else
            {
                lwork = g_lwork;
            }

            create_vector(datatype, &work, lwork);

            /* Call to ZGEHRD API to generate hessenberg matrix*/
            fla_lapack_zgehrd(&n, ilo, ihi, A, &lda, tau, work, &lwork, info);
            reset_vector(datatype, work, lwork, 1);
            copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

            /* Call to ZUNGHR API to generate orthogonal matrix*/
            fla_lapack_zunghr(&n, ilo, ihi, A_save, &lda, tau, work, &lwork, info);
            copy_matrix(datatype, "full", n, n, A_save, lda, Z, ldz);

            /* Convert matrix from ZGEHRD to upper hessenberg matrix */
            convert_upper_hessenberg(datatype, n, A, lda);

            free_vector(work);
            break;
        }
    }
    free_matrix(A_save);
    free_vector(tau);
}

/* Convert matrix to upper hessenberg form */
void convert_upper_hessenberg(integer datatype, integer n, void *A, integer lda)
{
    integer i;
    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < n; i++)
            {
                float *p = &((float *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < n; i++)
            {
                double *p = &((double *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < n; i++)
            {
                scomplex *p = &((scomplex *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < n; i++)
            {
                dcomplex *p = &((dcomplex *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
    }
}

/* Pack a symmetric matrix in column first order */
void pack_matrix_lt(integer datatype, void *A, void *B, integer N, integer lda)
{
    integer i, j;

    switch(datatype)
    {
        case FLOAT:
        {
            float *bptr = (float *)B;

            for(i = 0; i < N; i++)
            {
                for(j = i; j < N; j++)
                {
                    *bptr++ = ((float *)A)[i * lda + j];
                }
            }
            break;
        }
        case DOUBLE:
        {
            double *bptr = B;

            for(i = 0; i < N; i++)
            {
                for(j = i; j < N; j++)
                {
                    *bptr++ = ((double *)A)[i * lda + j];
                }
            }
            break;
        }
        case COMPLEX:
        {
            scomplex *bptr = B;

            for(i = 0; i < N; i++)
            {
                for(j = i; j < N; j++)
                {
                    bptr->real = ((scomplex *)A)[i * lda + j].real;
                    bptr->imag = ((scomplex *)A)[i * lda + j].imag;
                    bptr++;
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex *bptr = B;

            for(i = 0; i < N; i++)
            {
                for(j = i; j < N; j++)
                {
                    bptr->real = ((dcomplex *)A)[i * lda + j].real;
                    bptr->imag = ((dcomplex *)A)[i * lda + j].imag;
                    bptr++;
                }
            }
            break;
        }
    }
}

/* Convert matrix to upper hessenberg form */
void extract_upper_hessenberg_matrix(integer datatype, integer n, void *A, integer lda)
{
    integer i;
    switch(datatype)
    {
        case FLOAT:
        {
            /* Making elements below sub diagonal to Zero */
            for(i = 0; i < n; i++)
            {
                float *p = &((float *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
        case DOUBLE:
        {
            /* Making elements below sub diagonal to Zero */
            for(i = 0; i < n; i++)
            {
                double *p = &((double *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
        case COMPLEX:
        {
            /* Making elements below sub diagonal to Zero */
            for(i = 0; i < n; i++)
            {
                scomplex *p = &((scomplex *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Making elements below sub diagonal to Zero */
            for(i = 0; i < n; i++)
            {
                dcomplex *p = &((dcomplex *)A)[(i + 2) + i * lda];
                reset_vector(datatype, (void *)p, n - i - 2, 1);
            }
            break;
        }
    }
}

/* Decompose matrix A in to QR and store orthogonal matrix in Q and R in A*/
void get_orthogonal_matrix_from_QR(integer datatype, integer n, void *A, integer lda, void *Q,
                                   integer ldq, integer *info)
{
    void *tau = NULL, *work = NULL;
    integer lwork = -1;
    if((lda < n) || (ldq < n))
        return;
    /* Intializing matrix for the call to GGHRD */
    create_vector(datatype, &work, 1);
    create_vector(datatype, &tau, n);

    switch(datatype)
    {
        case FLOAT:
        {
            /* Generating orthogonal matrix Q by QR reduction of A */
            copy_matrix(datatype, "full", n, n, A, lda, Q, ldq);
            fla_lapack_sgeqrf(&n, &n, NULL, &ldq, NULL, work, &lwork, info);
            if(*info < 0)
                break;
            else
                lwork = get_work_value(datatype, work);
            free_vector(work);
            create_vector(datatype, &work, lwork);
            /* Call to SGEQRF to decompose matrix to QR form */
            fla_lapack_sgeqrf(&n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            reset_matrix(datatype, n, n, A, lda);
            copy_matrix(datatype, "Upper", n, n, Q, ldq, A, lda);
            reset_vector(datatype, work, lwork, 1);
            /* Call to SORGQR to calculate matrix to Q */
            fla_lapack_sorgqr(&n, &n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            break;
        }
        case DOUBLE:
        {
            /* Generating orthogonal matrix Q by QR reduction of A */
            copy_matrix(datatype, "full", n, n, A, lda, Q, ldq);
            fla_lapack_dgeqrf(&n, &n, NULL, &ldq, NULL, work, &lwork, info);
            if(*info < 0)
                break;
            else
                lwork = get_work_value(datatype, work);
            free_vector(work);
            create_vector(datatype, &work, lwork);
            /* Call to DGEQRF to decompose matrix to QR form */
            fla_lapack_dgeqrf(&n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            reset_matrix(datatype, n, n, A, lda);
            copy_matrix(datatype, "Upper", n, n, Q, ldq, A, lda);
            reset_vector(datatype, work, lwork, 1);
            /* Call to DORGQR to calculate matrix to Q */
            fla_lapack_dorgqr(&n, &n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            break;
        }
        case COMPLEX:
        {
            /* Generating orthogonal matrix Q by QR reduction of A */
            copy_matrix(datatype, "full", n, n, A, lda, Q, ldq);
            fla_lapack_cgeqrf(&n, &n, NULL, &ldq, NULL, work, &lwork, info);
            if(*info < 0)
                break;
            else
                lwork = get_work_value(datatype, work);
            free_vector(work);
            create_vector(datatype, &work, lwork);
            /* Call to CGEQRF to decompose matrix to QR form */
            fla_lapack_cgeqrf(&n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            reset_matrix(datatype, n, n, A, lda);
            copy_matrix(datatype, "Upper", n, n, Q, ldq, A, lda);
            reset_vector(datatype, work, lwork, 1);
            /* Call to CUNGQR to calculate matrix to Q */
            fla_lapack_cungqr(&n, &n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Generating orthogonal matrix Q by QR reduction of A */
            copy_matrix(datatype, "full", n, n, A, lda, Q, ldq);
            fla_lapack_zgeqrf(&n, &n, NULL, &ldq, NULL, work, &lwork, info);
            if(*info < 0)
                break;
            else
                lwork = get_work_value(datatype, work);
            free_vector(work);
            create_vector(datatype, &work, lwork);
            /* Call to ZGEQRF to decompose matrix to QR form */
            fla_lapack_zgeqrf(&n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            reset_matrix(datatype, n, n, A, lda);
            copy_matrix(datatype, "Upper", n, n, Q, ldq, A, lda);
            reset_vector(datatype, work, lwork, 1);
            /* Call to ZUNGQR to calculate matrix to Q */
            fla_lapack_zungqr(&n, &n, &n, Q, &ldq, tau, work, &lwork, info);
            if(*info < 0)
                break;
            break;
        }
    }
    free_vector(tau);
    free_vector(work);
}

/* Print matrix contents for visual inspection
 * if order == 'C' matrix will be printed in columns first order
 * else if order == 'R' matrix will be printed in rows first order
 */
void print_matrix(char *desc, char *order, integer datatype, integer M, integer N, void *A,
                  integer lda)
{
    integer i, j, row_max = M, col_max = N, ldc = lda, ldr = 1;
    if(*order == 'C')
    {
        row_max = N;
        col_max = M;
        ldc = 1;
        ldr = lda;
    }
    printf("\n %s:\n", desc);
    switch(datatype)
    {
        case INTEGER:
        {
            for(i = 0; i < row_max; i++)
            {
                for(j = 0; j < col_max; j++)
                {
                    printf(" %" FT_IS " ", ((integer *)A)[i * ldr + j * ldc]);
                }
                printf("\n");
            }
            break;
        }
        case FLOAT:
        {
            for(i = 0; i < row_max; i++)
            {
                for(j = 0; j < col_max; j++)
                {
                    printf(" %e", ((float *)A)[i * ldr + j * ldc]);
                }
                printf("\n");
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < row_max; i++)
            {
                for(j = 0; j < col_max; j++)
                {
                    printf(" %e", ((double *)A)[i * ldr + j * ldc]);
                }
                printf("\n");
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < row_max; i++)
            {
                for(j = 0; j < col_max; j++)
                {
                    printf(" (%e + j %e)", ((scomplex *)A)[i * ldr + j * ldc].real,
                           ((scomplex *)A)[i * ldr + j * ldc].imag);
                }
                printf("\n");
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < row_max; i++)
            {
                for(j = 0; j < col_max; j++)
                {
                    printf(" (%e + j %e)", ((dcomplex *)A)[i * ldr + j * ldc].real,
                           ((dcomplex *)A)[i * ldr + j * ldc].imag);
                }
                printf("\n");
            }
            break;
        }
    }
}

/* Get upper triangular matrix or lower triangular matrix based on UPLO */
void get_triangular_matrix(char *uplo, integer datatype, integer m, integer n, void *A, integer lda)
{
    rand_matrix(datatype, A, m, n, lda);
    integer i;

    switch(datatype)
    {
        case FLOAT:
        {
            if(*uplo == 'U')
            {
                for(i = 0; i < n; i++)
                {
                    float *p = &((float *)A)[(i + 1) + i * lda];
                    reset_vector(datatype, (void *)p, m - i - 1, 1);
                }
            }
            if(*uplo == 'L')
            {
                for(i = 0; i < n; i++)
                {
                    float *p = &((float *)A)[i * lda];
                    reset_vector(datatype, (void *)p, i + 1, 1);
                }
            }
            break;
        }
        case DOUBLE:
        {
            if(*uplo == 'U')
            {
                for(i = 0; i < n; i++)
                {
                    double *p = &((double *)A)[(i + 1) + i * lda];
                    reset_vector(datatype, (void *)p, m - i - 1, 1);
                }
            }
            if(*uplo == 'L')
            {
                for(i = 0; i < n; i++)
                {
                    double *p = &((double *)A)[i * lda];
                    reset_vector(datatype, (void *)p, i + 1, 1);
                }
            }
            break;
        }
        case COMPLEX:
        {
            if(*uplo == 'U')
            {
                for(i = 0; i < n; i++)
                {
                    scomplex *p = &((scomplex *)A)[(i + 1) + i * lda];
                    reset_vector(datatype, (void *)p, m - i - 1, 1);
                }
            }
            if(*uplo == 'L')
            {
                for(i = 0; i < n; i++)
                {
                    scomplex *p = &((scomplex *)A)[i * lda];
                    reset_vector(datatype, (void *)p, i + 1, 1);
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            if(*uplo == 'U')
            {
                for(i = 0; i < n; i++)
                {
                    dcomplex *p = &((dcomplex *)A)[(i + 1) + i * lda];
                    reset_vector(datatype, (void *)p, m - i - 1, 1);
                }
            }
            if(*uplo == 'L')
            {
                for(i = 0; i < n; i++)
                {
                    dcomplex *p = &((dcomplex *)A)[i * lda];
                    reset_vector(datatype, (void *)p, i + 1, 1);
                }
            }
            break;
        }
    }
}

/*Test to Check order of Singular values of SVD (positive and non-decreasing)*/
double svd_check_order(integer datatype, void *s, integer m, integer n, double residual)
{
    integer min_m_n, i;
    min_m_n = fla_min(m, n);
    double resid = 0.;

    switch(datatype)
    {
        case INTEGER:
        {
            for(i = 0; i < (min_m_n - 1); i++)
            {
                if((((int *)s)[i] < 0) || (((int *)s)[i] < ((int *)s)[i + 1]))
                {
                    resid = residual * 2;
                    break;
                }
            }
            if(((int *)s)[min_m_n - 1] < 0)
                resid = residual * 2;
            break;
        }
        case FLOAT:
        {
            for(i = 0; i < (min_m_n - 1); i++)
            {
                if((((float *)s)[i] < 0.f) || (((float *)s)[i] < ((float *)s)[i + 1]))
                {
                    resid = residual * 2;
                    break;
                }
            }
            if(((float *)s)[min_m_n - 1] < 0.f)
                resid = residual * 2;
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < (min_m_n - 1); i++)
            {
                if((((double *)s)[i] < 0.) || (((double *)s)[i] < ((double *)s)[i + 1]))
                {
                    resid = residual * 2;
                    break;
                }
            }
            if(((double *)s)[min_m_n - 1] < 0.)
                resid = residual * 2;
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < (min_m_n - 1); i++)
            {
                if((((float *)s)[i] < 0.f) || (((float *)s)[i] < ((float *)s)[i + 1]))
                {
                    resid = residual * 2;
                    break;
                }
            }
            if(((float *)s)[min_m_n - 1] < 0.f)
                resid = residual * 2;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < (min_m_n - 1); i++)
            {
                if((((double *)s)[i] < 0.) || (((double *)s)[i] < ((double *)s)[i + 1]))
                {
                    resid = residual * 2;
                    break;
                }
            }
            if(((double *)s)[min_m_n - 1] < 0.)
                resid = residual * 2;
            break;
        }
        default:
            break;
    }
    return resid;
}

/* Intialize matrix with special values*/
void init_matrix_spec_in(integer datatype, void *A, integer M, integer N, integer LDA, char type)
{
    integer i, j;
    if(LDA < M)
        return;
    switch(datatype)
    {
        case FLOAT:
        {
            float value = 0.f;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((float *)A)[i * LDA + j] = value;
                }
            }
            break;
        }
        case DOUBLE:
        {
            double value = 0.;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((double *)A)[i * LDA + j] = value;
                }
            }
            break;
        }
        case COMPLEX:
        {
            float value = 0.f;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((scomplex *)A)[i * LDA + j].real = value;
                    ((scomplex *)A)[i * LDA + j].imag = value;
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double value = 0.;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < N; i++)
            {
                for(j = 0; j < M; j++)
                {
                    ((dcomplex *)A)[i * LDA + j].real = value;
                    ((dcomplex *)A)[i * LDA + j].imag = value;
                }
            }
            break;
        }
    }
}

/* Intialize matrix with special values in random locations */
void init_matrix_spec_rand_in(integer datatype, void *A, integer M, integer N, integer LDA,
                              char type)
{
    integer rows, cols, upspan, lowspan, span;
    if(LDA < M)
        return;
    rand_matrix(datatype, A, M, N, LDA);
    /* when M*N less than 2 there is no need of randomness*/
    if(M * N < 2)
    {
        char type_;
        type_ = (type == 'A') ? 'N' : 'I';
        init_matrix_spec_in(datatype, A, M, N, LDA, type_);
        return;
    }
    /*
    Add random extreme values:
    for small size matrices, when M*N less than 10 adding one extreme value in upper triangular
    matrix, other one extreme value in lower triangular matrix
    for medium/large sizes, when M*N greater than 10 adding 10% of input values as extreme values
    in upper triangular matrix, other 10% of input values as exterme values in lower triangular
    matrix
    */
    if(M * N > 10)
    {
        lowspan = (M * N) * 0.1;
        upspan = (M * N) * 0.1;
    }
    else
    {
        lowspan = 1;
        upspan = 1;
    }
    span = lowspan + upspan;
    switch(datatype)
    {
        case FLOAT:
        {
            float value = 0.f;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                cols = rand() % N;
                /* Replace 10 percent of special values in upper triangular matrix */
                if(upspan > 0)
                {
                    if(rows <= cols)
                    {
                        if(!isnan(((float *)A)[cols * LDA + rows]))
                        {
                            ((float *)A)[cols * LDA + rows] = value;
                            upspan = upspan - 1;
                        }
                    }
                }
                /* Replace 10 percent of special values in lower triangular matrix */
                else if(lowspan > 0)
                {
                    if(rows >= cols)
                    {
                        if(!isnan(((float *)A)[cols * LDA + rows]))
                        {
                            ((float *)A)[cols * LDA + rows] = value;
                            lowspan = lowspan - 1;
                        }
                    }
                }
                span = lowspan + upspan;
            }
            break;
        }
        case DOUBLE:
        {
            double value = 0.;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                cols = rand() % N;
                if(upspan > 0)
                {
                    if(rows <= cols)
                    {
                        if(!isnan(((double *)A)[cols * LDA + rows]))
                        {
                            ((double *)A)[cols * LDA + rows] = value;
                            upspan = upspan - 1;
                        }
                    }
                }
                else if(lowspan > 0)
                {
                    if(rows >= cols)
                    {
                        if(!isnan(((double *)A)[cols * LDA + rows]))
                        {
                            ((double *)A)[cols * LDA + rows] = value;
                            lowspan = lowspan - 1;
                        }
                    }
                }
                span = lowspan + upspan;
            }
            break;
        }
        case COMPLEX:
        {
            float value = 0.f;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                cols = rand() % N;
                if(upspan > 0)
                {
                    if(rows <= cols)
                    {
                        if(!isnan(((scomplex *)A)[cols * LDA + rows].real))
                        {
                            ((scomplex *)A)[cols * LDA + rows].real = value;
                            ((scomplex *)A)[cols * LDA + rows].imag = value;
                            upspan = upspan - 1;
                        }
                    }
                }
                else if(lowspan > 0)
                {
                    if(rows >= cols)
                    {
                        if(!isnan(((scomplex *)A)[cols * LDA + rows].real))
                        {
                            ((scomplex *)A)[cols * LDA + rows].real = value;
                            ((scomplex *)A)[cols * LDA + rows].imag = value;
                            lowspan = lowspan - 1;
                        }
                    }
                }
                span = lowspan + upspan;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double value = 0.;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                cols = rand() % N;
                if(upspan > 0)
                {
                    if(rows <= cols)
                    {
                        if(!isnan(((dcomplex *)A)[cols * LDA + rows].real))
                        {
                            ((dcomplex *)A)[cols * LDA + rows].real = value;
                            ((dcomplex *)A)[cols * LDA + rows].imag = value;
                            upspan = upspan - 1;
                        }
                    }
                }
                else if(lowspan > 0)
                {
                    if(rows >= cols)
                    {
                        if(!isnan(((dcomplex *)A)[cols * LDA + rows].real))
                        {
                            ((dcomplex *)A)[cols * LDA + rows].real = value;
                            ((dcomplex *)A)[cols * LDA + rows].imag = value;
                            lowspan = lowspan - 1;
                        }
                    }
                }
                span = lowspan + upspan;
            }
            break;
        }
    }
}

/* Test to check the extreme values propagation in output matrix */
bool check_extreme_value(integer datatype, integer M, integer N, void *A, integer LDA, char type)
{
    if(!A)
        return false;
    integer i, j;
    switch(datatype)
    {
        case FLOAT:
        {
            if(type == 'A' || type == 'N')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if(isnan(((float *)A)[i * LDA + j]))
                        {
                            return true;
                        }
                    }
                }
            }

            else if(type == 'F' || type == 'I')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if((isinf(((float *)A)[i * LDA + j])) || (isnan(((float *)A)[i * LDA + j])))
                        {
                            return true;
                        }
                    }
                }
            }
            break;
        }

        case DOUBLE:
        {
            if(type == 'A' || type == 'N')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if(isnan(((double *)A)[i * LDA + j]))
                        {
                            return true;
                        }
                    }
                }
            }

            if(type == 'F' || type == 'I')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if((isinf(((double *)A)[i * LDA + j]))
                           || (isnan(((double *)A)[i * LDA + j])))
                        {
                            return true;
                        }
                    }
                }
            }
            break;
        }

        case COMPLEX:
        {
            if(type == 'A' || type == 'N')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if(isnan(((scomplex *)A)[i * LDA + j].real)
                           || isnan(((scomplex *)A)[i * LDA + j].imag))
                        {
                            return true;
                        }
                    }
                }
            }

            else if(type == 'F' || type == 'I')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if((isinf(((scomplex *)A)[i * LDA + j].real)
                            || isinf(((scomplex *)A)[i * LDA + j].imag))
                           || (isnan(((scomplex *)A)[i * LDA + j].real)
                               || isnan(((scomplex *)A)[i * LDA + j].imag)))
                        {
                            return true;
                        }
                    }
                }
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            if(type == 'A' || type == 'N')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if(isnan(((dcomplex *)A)[i * LDA + j].real)
                           || isnan(((dcomplex *)A)[i * LDA + j].imag))
                        {
                            return true;
                        }
                    }
                }
            }

            else if(type == 'F' || type == 'I')
            {
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < M; j++)
                    {
                        if((isinf(((dcomplex *)A)[i * LDA + j].real)
                            || isinf(((dcomplex *)A)[i * LDA + j].imag))
                           || (isnan(((dcomplex *)A)[i * LDA + j].real)
                               || isnan(((dcomplex *)A)[i * LDA + j].imag)))
                        {
                            return true;
                        }
                    }
                }
            }
            break;
        }
    }
    return 0;
}

/* Intialize vector with special values */
void init_vector_spec_in(integer datatype, void *A, integer M, integer incx, char type)
{
    integer i;
    switch(datatype)
    {
        case FLOAT:
        {
            float value = 0.f;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < M; i++)
            {
                ((float *)A)[i * incx] = value;
            }
            break;
        }
        case DOUBLE:
        {
            double value = 0.;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < M; i++)
            {
                ((double *)A)[i * incx] = value;
            }
            break;
        }
        case COMPLEX:
        {
            float value = 0.f;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < M; i++)
            {
                ((scomplex *)A)[i * incx].real = value;
                ((scomplex *)A)[i * incx].imag = value;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double value = 0.;
            if(type == 'I')
                value = INFINITY;
            else if(type == 'N')
                value = NAN;
            for(i = 0; i < M; i++)
            {
                ((dcomplex *)A)[i * incx].real = value;
                ((dcomplex *)A)[i * incx].imag = value;
            }
            break;
        }
    }
}
/* Intialize vector with special values in random locations */
void init_vector_spec_rand_in(integer datatype, void *A, integer M, integer incx, char type)
{
    integer rows, span;
    rand_vector(datatype, M, A, incx, d_zero, d_zero, 'R');
    /* when M*N less than 2 there is no need of randomness*/
    if(M < 2)
    {
        char type_;
        type_ = (type == 'A') ? 'N' : 'I';
        init_vector_spec_in(datatype, A, M, incx, type_);
        return;
    }
    /*
    Add random extreme values:
    for small size vector, when M  less than 10 adding two extreme value in vector,
    for medium/large sizes, when M greater than 10 adding 20% of input values as extreme values
    in vector
    */
    if(M > 10)
    {
        span = (M)*0.2;
    }
    else
    {
        span = 2;
    }
    switch(datatype)
    {
        case FLOAT:
        {
            float value = 0.f;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                /* Replace 20 percent of special values in vector */
                if(span > 0)
                {
                    if(!isnan(((float *)A)[rows * incx]))
                    {
                        ((float *)A)[rows * incx] = value;
                        span = span - 1;
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            double value = 0.;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                /* Replace 20 percent of special values in vector */
                if(span > 0)
                {
                    if(!isnan(((double *)A)[rows * incx]))
                    {
                        ((double *)A)[rows * incx] = value;
                        span = span - 1;
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            float value = 0.f;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                /* Replace 20 percent of special values in vector */
                if(span > 0)
                {
                    if(!isnan(((scomplex *)A)[rows * incx].real))
                    {
                        ((scomplex *)A)[rows * incx].real = value;
                        ((scomplex *)A)[rows * incx].imag = value;
                        span = span - 1;
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double value = 0.;
            if(type == 'F')
                value = INFINITY;
            else if(type == 'A')
                value = NAN;
            while(span > 0)
            {
                rows = rand() % M;
                /* Replace 20 percent of special values in vector */
                if(span > 0)
                {
                    if(!isnan(((dcomplex *)A)[rows * incx].real))
                    {
                        ((dcomplex *)A)[rows * incx].real = value;
                        ((dcomplex *)A)[rows * incx].imag = value;
                        span = span - 1;
                    }
                }
            }
            break;
        }
    }

    return;
}

/*Intialize matrix according to given input*/
void init_matrix(integer datatype, void *A, integer M, integer N, integer LDA, FILE *g_ext_fptr,
                 char imatrix_char)
{
    if(g_ext_fptr != NULL)
        init_matrix_from_file(datatype, A, M, N, LDA, g_ext_fptr);
    else if(imatrix_char == 'I' || imatrix_char == 'N')
        init_matrix_spec_in(datatype, A, M, N, LDA, imatrix_char);
    else if(imatrix_char == 'A' || imatrix_char == 'F')
        init_matrix_spec_rand_in(datatype, A, M, N, LDA, imatrix_char);
    else
        rand_matrix(datatype, A, M, N, LDA);
}

/*
 *   Create input matrix A by randomly generating singular values(S)
 *   according to range A, I, V respectively
 *   A -> All singular values => rand vector
 *   V -> Singular value between (vl, vu) range => vector in range(vl, vu)
 *   I -> Singular value between (il, iu) index
 *   U -> Singular value between (vl, vu) range with uniformly distributed values
 *               A  = (U * S * V')
 *   where  S  is a diagonal matrix with diagonal elements being
 *   the singular values of input matrix A
 *   U is an M-by-M orthogonal matrix, and
 *   V is an N-by-N orthogonal matrix.
 */
void create_svd_matrix(integer datatype, char range, integer m, integer n, void *A_input,
                       integer lda, void *S, double vl, double vu, integer il, integer iu,
                       char imatrix, void *scal, integer info)
{
    if(lda < m)
        return;
    /* For range I, index range should be 1 <= IL <= IU <= min(M,N) */
    if((range == 'I' || range == 'i') && (il <= 0 || iu > fla_min(m, n)))
        return;

    void *A, *B, *U, *V, *sigma, *Usigma, *s_test;
    integer min_m_n = fla_min(m, n);

    create_matrix(datatype, &A, m, m);
    create_matrix(datatype, &B, n, n);
    /* Orthogonal matrix U and V */
    create_matrix(datatype, &U, m, m);
    create_matrix(datatype, &V, n, n);
    /* Singular values array s_test */
    create_realtype_vector(datatype, &s_test, min_m_n);

    /* Generate random matrix for Decomposition */
    rand_matrix(datatype, A, m, m, m);
    rand_matrix(datatype, B, n, n, n);

    /* Calculate orthogonal matrix U & V */
    get_orthogonal_matrix_from_QR(datatype, m, A, m, U, m, &info);
    get_orthogonal_matrix_from_QR(datatype, n, B, n, V, n, &info);

    /* Generating positive singular values according to the ranges */
    rand_vector(get_realtype(datatype), min_m_n, s_test, i_one, vl, vu, range);

    /* Sorting singular values in descending order */
    get_abs_vector_value(datatype, s_test, min_m_n, i_one);
    sort_realtype_vector(datatype, "D", min_m_n, s_test, i_one);

    /* Copying the singular value to S array with respect to range for validation */
    if(range == 'I' || range == 'i')
    {
        copy_realtype_subvector(datatype, (iu - il + 1), s_test, S, (il - 1));
    }
    else
    {
        copy_realtype_vector(datatype, min_m_n, s_test, 1, S, 1);
    }

    /* Generating A matrix by A = (U * Sigma * VT) */
    create_matrix(datatype, &Usigma, m, n);
    create_matrix(datatype, &sigma, m, n);

    /* Diagonalize the singular values for sigma */
    diagonalize_realtype_vector(datatype, s_test, sigma, min_m_n, n, min_m_n);

    reset_matrix(datatype, m, n, Usigma, m);
    reset_matrix(datatype, m, n, A_input, lda);

    switch(datatype)
    {
        case FLOAT:
        {
            sgemm_("N", "N", &m, &n, &min_m_n, &s_one, U, &m, sigma, &min_m_n, &s_zero, Usigma, &m);
            sgemm_("N", "N", &m, &n, &n, &s_one, Usigma, &m, V, &n, &s_zero, A_input, &lda);
            break;
        }
        case DOUBLE:
        {
            dgemm_("N", "N", &m, &n, &min_m_n, &d_one, U, &m, sigma, &min_m_n, &d_zero, Usigma, &m);
            dgemm_("N", "N", &m, &n, &n, &d_one, Usigma, &m, V, &n, &d_zero, A_input, &lda);
            break;
        }
        case COMPLEX:
        {
            cgemm_("N", "N", &m, &n, &min_m_n, &c_one, U, &m, sigma, &min_m_n, &c_zero, Usigma, &m);
            cgemm_("N", "N", &m, &n, &n, &c_one, Usigma, &m, V, &n, &c_zero, A_input, &lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zgemm_("N", "N", &m, &n, &min_m_n, &z_one, U, &m, sigma, &min_m_n, &z_zero, Usigma, &m);
            zgemm_("N", "N", &m, &n, &n, &z_one, Usigma, &m, V, &n, &z_zero, A_input, &lda);
            break;
        }
    }

    if(imatrix == 'O' || imatrix == 'U')
    {
        /* Initializing matrix with values around overflow underflow */
        init_matrix_overflow_underflow_svd(datatype, m, n, A_input, lda, imatrix, scal);
    }

    free_matrix(A);
    free_matrix(B);
    free_matrix(U);
    free_matrix(V);
    free_vector(s_test);
    free_matrix(sigma);
    free_matrix(Usigma);
}

/* Copying vector between specified ranges
 * index - Starting position of the source vector */
void copy_realtype_subvector(integer datatype, integer m, void *A, void *B, integer index)
{
    if(datatype == FLOAT || datatype == COMPLEX)
    {
        float *float_A, *float_B;
        float_A = (float *)A + index;
        float_B = (float *)B;
        copy_vector(datatype, m, float_A, 1, float_B, 1);
    }
    else if(datatype == DOUBLE || datatype == DOUBLE_COMPLEX)
    {
        double *double_A, *double_B;
        double_A = (double *)A + index;
        double_B = (double *)B;
        copy_vector(datatype, m, double_A, 1, double_B, 1);
    }
}

/* Checks whether the value is zero or not */
double is_value_zero(integer datatype, void *value, double residual)
{
    double resid = residual;
    switch(datatype)
    {
        case FLOAT:
        {
            if(*(float *)value != s_zero)
            {
                resid = residual * 2.0;
            }
            break;
        }
        case DOUBLE:
        {
            if(*(double *)value != d_zero)
            {
                resid = residual * 2.0;
            }
            break;
        }
        case COMPLEX:
        {
            if(((scomplex *)value)[0].real != s_zero || ((scomplex *)value)[0].imag != s_zero)
            {
                resid = residual * 2.0;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            if(((dcomplex *)value)[0].real != d_zero || ((dcomplex *)value)[0].imag != d_zero)
            {
                resid = residual * 2.0;
            }
            break;
        }
    }
    return resid;
}

/*Intialize vector according to given input*/
void init_vector(integer datatype, void *A, integer M, integer incx, FILE *g_ext_fptr,
                 char ivector_char)
{
    if(g_ext_fptr != NULL)
        init_vector_from_file(datatype, A, M, incx, g_ext_fptr);
    else if(ivector_char == 'I' || ivector_char == 'N')
        init_vector_spec_in(datatype, A, M, incx, ivector_char);
    else if(ivector_char == 'A' || ivector_char == 'F')
        init_vector_spec_rand_in(datatype, A, M, incx, ivector_char);
    else
        rand_vector(datatype, M, A, incx, d_zero, d_zero, 'R');
}

/* General matrix multiplication of tridiagonal matrix using LAGTM
 * C = tridiag_matrix[du, d, dl] * B */
void tridiag_matrix_multiply(integer datatype, integer n, integer nrhs, void *dl, void *d, void *du,
                             void *B, integer ldb, void *C, integer ldc)
{
    if(ldb < n || ldc < n)
        return;

    switch(datatype)
    {
        case FLOAT:
        {
            slagtm_("N", &n, &nrhs, &s_one, dl, d, du, B, &ldb, &s_zero, C, &ldc);
            break;
        }
        case DOUBLE:
        {
            dlagtm_("N", &n, &nrhs, &d_one, dl, d, du, B, &ldb, &d_zero, C, &ldc);
            break;
        }
        case COMPLEX:
        {
            clagtm_("N", &n, &nrhs, &s_one, dl, d, du, B, &ldb, &s_zero, C, &ldc);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zlagtm_("N", &n, &nrhs, &d_one, dl, d, du, B, &ldb, &d_zero, C, &ldc);
            break;
        }
    }
}

/* Calculate the difference between two matrix  A = A - B */
void matrix_difference(integer datatype, integer m, integer n, void *A, integer lda, void *B,
                       integer ldb)
{
    integer i = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            float *float_A, *float_B;
            for(i = 0; i < n; i++)
            {
                float_A = ((float *)A + (i * lda));
                float_B = ((float *)B + (i * ldb));
                saxpy_(&m, &s_n_one, float_B, &i_one, float_A, &i_one);
            }
            break;
        }
        case DOUBLE:
        {
            double *double_A, *double_B;
            for(i = 0; i < n; i++)
            {
                double_A = ((double *)A + (i * lda));
                double_B = ((double *)B + (i * ldb));
                daxpy_(&m, &d_n_one, double_B, &i_one, double_A, &i_one);
            }
            break;
        }
        case COMPLEX:
        {
            scomplex *scomplex_A, *scomplex_B;
            for(i = 0; i < n; i++)
            {
                scomplex_A = ((scomplex *)A + (i * lda));
                scomplex_B = ((scomplex *)B + (i * ldb));
                caxpy_(&m, &c_n_one, scomplex_B, &i_one, scomplex_A, &i_one);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex *dcomplex_A, *dcomplex_B;
            for(i = 0; i < n; i++)
            {
                dcomplex_A = ((dcomplex *)A + (i * lda));
                dcomplex_B = ((dcomplex *)B + (i * ldb));
                zaxpy_(&m, &z_n_one, dcomplex_B, &i_one, dcomplex_A, &i_one);
            }
            break;
        }
    }
}

/* Copy tridiagonal matrix */
void copy_tridiag_matrix(integer datatype, void *dl, void *d, void *du, integer M, integer N,
                         void *A, integer LDA)
{
    if(LDA < M)
        return;
    integer inc = LDA + 1;
    integer min_m_n, dl_size, du_size;
    min_m_n = fla_min(M, N);
    dl_size = N - 1, du_size = N - 1;

    reset_matrix(datatype, M, N, A, LDA);

    switch(datatype)
    {
        case FLOAT:
        {
            float *d_ptr = ((float *)A);
            float *du_ptr = ((float *)A + LDA);
            float *dl_ptr = ((float *)A + 1);
            scopy_(&min_m_n, d, &i_one, d_ptr, &inc);
            scopy_(&du_size, du, &i_one, du_ptr, &inc);
            scopy_(&dl_size, dl, &i_one, dl_ptr, &inc);
            break;
        }
        case DOUBLE:
        {
            double *d_ptr = ((double *)A);
            double *du_ptr = ((double *)A + LDA);
            double *dl_ptr = ((double *)A + 1);
            dcopy_(&min_m_n, d, &i_one, d_ptr, &inc);
            dcopy_(&du_size, du, &i_one, du_ptr, &inc);
            dcopy_(&dl_size, dl, &i_one, dl_ptr, &inc);
            break;
        }
        case COMPLEX:
        {
            scomplex *d_ptr = ((scomplex *)A);
            scomplex *du_ptr = ((scomplex *)A + LDA);
            scomplex *dl_ptr = ((scomplex *)A + 1);
            ccopy_(&min_m_n, d, &i_one, d_ptr, &inc);
            ccopy_(&du_size, du, &i_one, du_ptr, &inc);
            ccopy_(&dl_size, dl, &i_one, dl_ptr, &inc);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex *d_ptr = ((dcomplex *)A);
            dcomplex *du_ptr = ((dcomplex *)A + LDA);
            dcomplex *dl_ptr = ((dcomplex *)A + 1);
            zcopy_(&min_m_n, d, &i_one, d_ptr, &inc);
            zcopy_(&du_size, du, &i_one, du_ptr, &inc);
            zcopy_(&dl_size, dl, &i_one, dl_ptr, &inc);
            break;
        }
    }
}

/* Multiply general m * n matrix with diagonal real type vector
   (of an n * n diagonal matrix) of size n
   NOTE: General matrix by vector multiplication can be done by scaling each
         column of the matrix with corresponding element in the vector */
void multiply_matrix_diag_vector(integer datatype, integer m, integer n, void *A, integer lda,
                                 void *X, integer incx)
{
    integer j;
    if(m <= 0 || n <= 0)
        return;

    switch(datatype)
    {
        case FLOAT:
        {
            float *a_begin, *x = (float *)X;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by corresponding element
                   in the vector */
                a_begin = (float *)A + j * lda;
                sscal_(&m, &x[j], a_begin, &incx);
            }
            break;
        }

        case DOUBLE:
        {
            double *a_begin, *x = (double *)X;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by corresponding element
                   in the vector */
                a_begin = (double *)A + j * lda;
                dscal_(&m, &x[j], a_begin, &incx);
            }
            break;
        }

        case COMPLEX:
        {
            scomplex *a_begin, x;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by corresponding element
                   in the vector */
                a_begin = (scomplex *)A + j * lda;
                x.real = *((float *)X + j);
                x.imag = 0.f;
                cscal_(&m, &x, a_begin, &incx);
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            dcomplex *a_begin, x;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by corresponding element
                   in the vector */
                a_begin = (dcomplex *)A + j * lda;
                x.real = *((double *)X + j);
                x.imag = 0.0;
                zscal_(&m, &x, a_begin, &incx);
            }
            break;
        }
    }
}

/*
 * Generate square matrix of size n x n using Eigen decomposition(ED)
 *                     A  = (Q * lambda * Q')
 * where Q is an n x n orthogonal matrix and
 *       lambda is a diagonal vector(realtype)
 * NOTE: For simplification of general matrix - diagonal matrix multiplication,
 *       in this funciton lamda is taken as a vector.
 */
void generate_matrix_from_ED(integer datatype, integer n, void *A, integer lda, void *Q,
                             void *lambda)
{
    void *Qlambda = NULL;
    if(lda < n)
        return;

    create_matrix(datatype, &Qlambda, n, n);
    copy_matrix(datatype, "full", n, n, Q, n, Qlambda, n);

    /* Perform Q * lambda */
    multiply_matrix_diag_vector(datatype, n, n, Qlambda, n, lambda, 1);
    switch(datatype)
    {
        case FLOAT:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            sgemm_("N", "T", &n, &n, &n, &s_one, Qlambda, &n, Q, &n, &s_zero, A, &lda);
            break;
        }
        case DOUBLE:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            dgemm_("N", "T", &n, &n, &n, &d_one, Qlambda, &n, Q, &n, &d_zero, A, &lda);
            break;
        }
        case COMPLEX:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            cgemm_("N", "C", &n, &n, &n, &c_one, Qlambda, &n, Q, &n, &c_zero, A, &lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            zgemm_("N", "C", &n, &n, &n, &z_one, Qlambda, &n, Q, &n, &z_zero, A, &lda);
            break;
        }
    }
    free_matrix(Qlambda);
}

/* Sort the given real type vector in the given order
   order = A - Ascending order
         = D - Descending order */
void sort_realtype_vector(integer datatype, char *order, integer vect_len, void *w, integer incw)
{
    if(!w)
        return;

    integer i, j;
    if(get_realtype(datatype) == FLOAT)
    {
        float temp;
        float *w_ptr = (float *)w;

        for(i = 0; i < vect_len; i++)
        {
            for(j = i + 1; j < vect_len; j++)
            {
                if(*order == 'A')
                {
                    if(*(w_ptr + i * incw) > *(w_ptr + j * incw))
                    {
                        temp = *(w_ptr + i * incw);
                        *(w_ptr + i * incw) = *(w_ptr + j * incw);
                        *(w_ptr + j * incw) = temp;
                    }
                }
                else if(*order == 'D')
                {
                    if(*(w_ptr + i * incw) < *(w_ptr + j * incw))
                    {
                        temp = *(w_ptr + i * incw);
                        *(w_ptr + i * incw) = *(w_ptr + j * incw);
                        *(w_ptr + j * incw) = temp;
                    }
                }
            }
        }
    }
    else if(get_realtype(datatype) == DOUBLE)
    {
        double temp;
        double *w_ptr = (double *)w;

        for(i = 0; i < vect_len; i++)
        {
            for(j = i + 1; j < vect_len; j++)
            {
                if(*order == 'A')
                {
                    if(*(w_ptr + i * incw) > *(w_ptr + j * incw))
                    {
                        temp = *(w_ptr + i * incw);
                        *(w_ptr + i * incw) = *(w_ptr + j * incw);
                        *(w_ptr + j * incw) = temp;
                    }
                }
                else if(*order == 'D')
                {
                    if(*(w_ptr + i * incw) < *(w_ptr + j * incw))
                    {
                        temp = *(w_ptr + i * incw);
                        *(w_ptr + i * incw) = *(w_ptr + j * incw);
                        *(w_ptr + j * incw) = temp;
                    }
                }
            }
        }
    }
}

/* Compare two vectors starting from offset_A in A vector with B vector
   (starting from offset 0 in B) */
integer compare_realtype_vector(integer datatype, integer vect_len, void *A, integer inca,
                                integer offset_A, void *B, integer incb)
{
    if(!A || !B)
        return 1;

    integer i;
    if(datatype == FLOAT || datatype == COMPLEX)
    {
        float *a = (float *)A;
        float *b = (float *)B;
        for(i = 0; i < vect_len; i++)
        {
            if(f2c_abs(a[(i * inca) + (offset_A - 1)] - b[i * incb]) > MAX_FLT_DIFF)
            {
                return 1;
            }
        }
    }
    else
    {
        double *a = (double *)A;
        double *b = (double *)B;
        for(i = 0; i < vect_len; i++)
        {
            if(f2c_abs(a[(i * inca) + (offset_A - 1)] - b[i * incb]) > MAX_DBL_DIFF)
                return 1;
        }
    }
    return 0;
}

/*
 *   Create input matrix A by randomly generating eigen values(EVs) in given
 *   range (vl,vu)
 *               A  = (Q * lambda * Q')
 *   where  lambda  is a diagonal matrix with diagonal elements being
 *                  the eigen values of input matrix A
 *               Q  is an orthogonal matrix with corresponding
 *                  eigen vectors as its rows.
 */
void generate_matrix_from_EVs(integer datatype, char range, integer n, void *A, integer lda,
                              void *L, double vl, double vu)
{
    void *X = NULL, *Q = NULL;
    integer realtype, info = 0;
    realtype = get_realtype(datatype);

    /* Generate random vector of size n with values ranging
    between (vl, vu) */
    rand_vector(realtype, n, L, i_one, vl, vu, range);

    create_matrix(datatype, &X, n, n);
    rand_matrix(datatype, X, n, n, n);
    create_matrix(datatype, &Q, n, n);
    /* Generate random orthogonal matrix(Q) of size n x n */
    get_orthogonal_matrix_from_QR(datatype, n, X, n, Q, n, &info);
    /* Generate input matrix A using L(Eigen values)
       and Q(Eigen vectors) obtained above using reverse
       Eigen decompostion */
    generate_matrix_from_ED(datatype, n, A, lda, Q, L);
    /* Free up the buffers */
    free_matrix(X);
    free_matrix(Q);
}

/* Get absolute value of a real vector*/
void get_abs_vector_value(integer datatype, void *S, integer M, integer inc)
{
    integer i = 0;
    if(datatype == FLOAT || datatype == COMPLEX)
    {
        for(i = 0; i < M; i++)
            ((float *)S)[i * inc] = FLA_FABS(((float *)S)[i * inc]);
    }
    else if(datatype == DOUBLE || datatype == DOUBLE_COMPLEX)
    {
        for(i = 0; i < M; i++)
            ((double *)S)[i * inc] = FLA_FABS(((double *)S)[i * inc]);
    }
}

/* Initialize band matrix with random values.
Note: Input buffer A has to be allocated by caller.*/
void rand_band_matrix(integer datatype, integer M, integer N, integer kl, integer ku, void *A,
                      integer LDA)
{
    integer i, j, min_m_n;

    if((M <= 0) || (N <= 0) || (kl < 0) || (ku < 0) || (LDA <= 0) || (LDA < M) || (A == NULL))
        return;

    min_m_n = fla_min(M, N);
    reset_matrix(datatype, M, N, A, LDA);

    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < min_m_n; i++)
            {
                for(j = fla_max(0, i - kl); j < fla_min(min_m_n, i + ku + 1); j++)
                {
                    ((float *)A)[i + j * LDA] = SRAND();
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < min_m_n; i++)
            {
                for(j = fla_max(0, i - kl); j < fla_min(min_m_n, i + ku + 1); j++)
                {
                    ((double *)A)[i + j * LDA] = DRAND();
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < min_m_n; i++)
            {
                for(j = fla_max(0, i - kl); j < fla_min(min_m_n, i + ku + 1); j++)
                {
                    ((scomplex *)A)[i + j * LDA].real = SRAND();
                    ((scomplex *)A)[i + j * LDA].imag = SRAND();
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < min_m_n; i++)
            {
                for(j = fla_max(0, i - kl); j < fla_min(min_m_n, i + ku + 1); j++)
                {
                    ((dcomplex *)A)[i + j * LDA].real = DRAND();
                    ((dcomplex *)A)[i + j * LDA].imag = DRAND();
                }
            }
            break;
        }
    }
}

/* Initialize band storage for given band matrix.
Note: Input buffer A has to be allocated, initialized with band matrix by caller.*/
void get_band_storage_matrix(integer datatype, integer M, integer N, integer kl, integer ku,
                             void *A, integer LDA, void *AB, integer LDAB)
{
    integer i, j, A_size = LDA * N, AB_size = LDAB * N;

    if((M <= 0) || (N <= 0) || (kl < 0) || (ku < 0) || (LDA <= 0) || (LDA < M) || (A == NULL)
       || (AB == NULL) || (LDAB <= 0) || (LDAB < (2 * kl + ku + 1)))
        return;

    reset_matrix(datatype, LDAB, N, AB, LDAB);

    switch(datatype)
    {
        case FLOAT:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((float *)AB)[(kl + ku + i - j) + j * LDAB] = ((float *)A)[i + j * LDA];
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((double *)AB)[(kl + ku + i - j) + j * LDAB] = ((double *)A)[i + j * LDA];
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((scomplex *)AB)[(kl + ku + i - j) + j * LDAB].real
                            = ((scomplex *)A)[i + j * LDA].real;
                        ((scomplex *)AB)[(kl + ku + i - j) + j * LDAB].imag
                            = ((scomplex *)A)[i + j * LDA].imag;
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((dcomplex *)AB)[(kl + ku + i - j) + j * LDAB].real
                            = ((dcomplex *)A)[i + j * LDA].real;
                        ((dcomplex *)AB)[(kl + ku + i - j) + j * LDAB].imag
                            = ((dcomplex *)A)[i + j * LDA].imag;
                    }
                }
            }
            break;
        }
    }
}

/* Get band matrix from band storage matrix.
Note: Input buffer A has to be allocated, initialized with band storage matrix by caller.*/
void get_band_matrix_from_band_storage(integer datatype, integer M, integer N, integer kl,
                                       integer ku, void *AB, integer LDAB, void *A, integer LDA)
{
    integer i, j, A_size = LDA * N, AB_size = LDAB * N;

    if((M <= 0) || (N <= 0) || (kl < 0) || (ku < 0) || (LDA <= 0) || (LDA < M) || (A == NULL)
       || (AB == NULL) || (LDAB <= 0) || (LDAB < (2 * kl + ku + 1)))
        return;
    reset_matrix(datatype, M, N, A, LDA);
    switch(datatype)
    {
        case FLOAT:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - kl - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((float *)A)[i + j * LDA] = ((float *)AB)[(kl + ku + i - j) + j * LDAB];
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((double *)A)[i + j * LDA] = ((double *)AB)[(kl + ku + i - j) + j * LDAB];
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((scomplex *)A)[i + j * LDA].real
                            = ((scomplex *)AB)[(kl + ku + i - j) + j * LDAB].real;
                        ((scomplex *)A)[i + j * LDA].imag
                            = ((scomplex *)AB)[(kl + ku + i - j) + j * LDAB].imag;
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(j = 0; j < N; j++)
            {
                for(i = fla_max(0, j - ku); i <= fla_min(M - 1, j + kl); i++)
                {
                    if((((kl + ku + i - j) + j * LDAB) < AB_size) && ((i + j * LDA) < A_size))
                    {
                        ((dcomplex *)A)[i + j * LDA].real
                            = ((dcomplex *)AB)[(kl + ku + i - j) + j * LDAB].real;
                        ((dcomplex *)A)[i + j * LDA].imag
                            = ((dcomplex *)AB)[(kl + ku + i - j) + j * LDAB].imag;
                    }
                }
            }
            break;
        }
    }
}

/* Initialize band storage with random band matrix.
   Note: Input buffer AB has to be allocated by caller.*/
void rand_band_storage_matrix(integer datatype, integer M, integer N, integer kl, integer ku,
                              void *AB, integer LDAB)
{
    void *A;

    /* Allocate matrix A */
    create_matrix(datatype, &A, M, N);

    /* Initialize rand band matrix */
    rand_band_matrix(datatype, M, N, kl, ku, A, M);

    /* Convert band matrix into band storage. */
    get_band_storage_matrix(datatype, M, N, kl, ku, A, M, AB, LDAB);

    free_matrix(A);
}

/* On input, AB is the output of GBTRF().
   On output, AB is the reconstructed band storage matrix same as input of GBTRF().*/
void reconstruct_band_storage_matrix(integer datatype, integer m, integer n, integer kl, integer ku,
                                     void *AB, integer ldab, integer *ipiv)
{
    void *ABfac;
    integer j, i__1, i__2, i;
    integer diag_offset, superdiag_band, column_length;
    integer il, ip, iw;

    /* Copy factorized banded storage matrix */
    create_matrix(datatype, &ABfac, ldab, n);
    copy_matrix(datatype, "full", ldab, n, AB, ldab, ABfac, ldab);

    /* Reset matrix A to store reconstructed band matrix*/
    reset_matrix(datatype, ldab, n, AB, ldab);

    /* Iterate over each column */
    diag_offset = kl + ku;
    for(j = 0; j < n; ++j)
    {
        /* Determine super and sub diagonal band size */
        superdiag_band = fla_min(kl + ku, j);

        /* Determine column length of upper diagonal band */
        column_length = fla_min(m - 1, j) - j + superdiag_band + 1;

        if(column_length > 0)
        {
            switch(datatype)
            {
                case FLOAT:
                {
                    float t;

                    /* Copy current column of upper diagonal band to corresponding column of matrix
                     * A*/
                    copy_matrix(datatype, "full", column_length, 1,
                                &((float *)ABfac)[diag_offset - superdiag_band + j * ldab], ldab,
                                &((float *)AB)[diag_offset - column_length + 1 + j * ldab], ldab);

                    /* Apply multipliers of lower diagonal band to compute sub diagonal band
                     * elements  */
                    i__1 = m - 1;
                    i__2 = j - superdiag_band;
                    for(i = fla_min(i__1, j); i >= i__2; --i)
                    {
                        il = fla_min(kl, m - i - 1);
                        if(il > 0)
                        {
                            iw = diag_offset - 1 + i + 1 - j + (j)*ldab;
                            t = ((float *)AB)[iw];
                            saxpy_(&il, &t, &((float *)ABfac)[diag_offset + 1 + i * ldab], &i_one,
                                   &((float *)AB)[iw + 1], &i_one);

                            /* Swap the elements of the current column with the pivot */
                            ip = ipiv[i] - 1;
                            if(i != ip)
                            {
                                ip = ip - j + superdiag_band
                                     + fla_max(diag_offset - superdiag_band, 0) + j * ldab;
                                ((float *)AB)[iw] = ((float *)AB)[ip];
                                ((float *)AB)[ip] = t;
                            }
                        }
                    }
                    break;
                }
                case DOUBLE:
                {
                    double t;
                    /* Copy current column of upper diagonal band to corresponding column of matrix
                     * A*/
                    copy_matrix(datatype, "full", column_length, 1,
                                &((double *)ABfac)[diag_offset - superdiag_band + j * ldab], ldab,
                                &((double *)AB)[diag_offset - column_length + 1 + j * ldab], ldab);

                    /* Apply multipliers of lower diagonal band to compute sub diagonal band
                     * elements  */
                    i__1 = m - 1;
                    i__2 = j - superdiag_band;
                    for(i = fla_min(i__1, j); i >= i__2; --i)
                    {
                        il = fla_min(kl, m - i - 1);
                        if(il > 0)
                        {
                            iw = diag_offset - 1 + i + 1 - j + (j)*ldab;
                            t = ((double *)AB)[iw];
                            daxpy_(&il, &t, &((double *)ABfac)[diag_offset + 1 + i * ldab], &i_one,
                                   &((double *)AB)[iw + 1], &i_one);

                            /* Swap the elements of the current column with the pivot */
                            ip = ipiv[i] - 1;
                            if(i != ip)
                            {
                                ip = ip - j + superdiag_band
                                     + fla_max(diag_offset - superdiag_band, 0) + j * ldab;
                                ((double *)AB)[iw] = ((double *)AB)[ip];
                                ((double *)AB)[ip] = t;
                            }
                        }
                    }
                    break;
                }
                case COMPLEX:
                {
                    scomplex t;

                    /* Copy current column of upper diagonal band to corresponding column of matrix
                     * A*/
                    copy_matrix(datatype, "full", column_length, 1,
                                &((scomplex *)ABfac)[diag_offset - superdiag_band + j * ldab], ldab,
                                &((scomplex *)AB)[diag_offset - column_length + 1 + j * ldab],
                                ldab);

                    /* Apply multipliers of lower diagonal band to compute sub diagonal band
                     * elements  */
                    i__1 = m - 1;
                    i__2 = j - superdiag_band;
                    for(i = fla_min(i__1, j); i >= i__2; --i)
                    {
                        il = fla_min(kl, m - i - 1);
                        if(il > 0)
                        {
                            iw = diag_offset - 1 + i + 1 - j + (j)*ldab;
                            t = ((scomplex *)AB)[iw];
                            caxpy_(&il, &t, &((scomplex *)ABfac)[diag_offset + 1 + i * ldab],
                                   &i_one, &((scomplex *)AB)[iw + 1], &i_one);

                            /* Swap the elements of the current column with the pivot */
                            ip = ipiv[i] - 1;
                            if(i != ip)
                            {
                                ip = ip - j + superdiag_band
                                     + fla_max(diag_offset - superdiag_band, 0) + j * ldab;
                                ((scomplex *)AB)[iw] = ((scomplex *)AB)[ip];
                                ((scomplex *)AB)[ip] = t;
                            }
                        }
                    }
                    break;
                }
                case DOUBLE_COMPLEX:
                {
                    dcomplex t;

                    /* Copy current column of upper diagonal band to corresponding column of matrix
                     * A*/
                    copy_matrix(datatype, "full", column_length, 1,
                                &((dcomplex *)ABfac)[diag_offset - superdiag_band + j * ldab], ldab,
                                &((dcomplex *)AB)[diag_offset - column_length + 1 + j * ldab],
                                ldab);

                    /* Apply multipliers of lower diagonal band to compute sub diagonal band
                     * elements  */
                    i__1 = m - 1;
                    i__2 = j - superdiag_band;
                    for(i = fla_min(i__1, j); i >= i__2; --i)
                    {
                        il = fla_min(kl, m - i - 1);
                        if(il > 0)
                        {
                            iw = diag_offset - 1 + i + 1 - j + (j)*ldab;
                            t = ((dcomplex *)AB)[iw];
                            zaxpy_(&il, &t, &((dcomplex *)ABfac)[diag_offset + 1 + i * ldab],
                                   &i_one, &((dcomplex *)AB)[iw + 1], &i_one);

                            /* Swap the elements of the current column with the pivot */
                            ip = ipiv[i] - 1;
                            if(i != ip)
                            {
                                ip = ip - j + superdiag_band
                                     + fla_max(diag_offset - superdiag_band, 0) + j * ldab;
                                ((dcomplex *)AB)[iw] = ((dcomplex *)AB)[ip];
                                ((dcomplex *)AB)[ip] = t;
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    free_matrix(ABfac);
}

/* Test for checking whether solution x of Ax = B from least square api belongs to
 * row space of A
 */
void check_vector_in_rowspace(integer datatype, char *trans, integer m, integer n, integer nrhs,
                              void *A, integer lda, void *x, integer ldb, void *resid)
{
    integer lwork = (m + nrhs) * (n + 2), ldwork = (m + nrhs), temp, INFO = 0;
    double temp1;
    void *work = NULL;

    if(*trans == 'T' || *trans == 'C')
    {
        lwork = (n + nrhs) * (m + 2);
        ldwork = m;
    }
    create_vector(datatype, &work, lwork);

    /* checks whether X is in the row space of A or A'.  It does so
     * by scaling both X and A such that their norms are in the range
     * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
     * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
     * and returning the norm of the trailing triangle, scaled by
     * MAX(M,N,NRHS)*eps.
     */
    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm = 0, norm_a = 0, norm_x = 0, rwork;
            eps = fla_lapack_slamch("E");

            if((*trans == 'T' && m > n) || (*trans == 'N' && m < n))
            {
                /* Copy A into work */
                fla_lapack_slacpy("All", &m, &n, A, &lda, work, &ldwork);
                norm_a = fla_lapack_slange("M", &m, &n, work, &ldwork, &rwork);
                /*Scale work*/
                if(norm_a != 0.)
                {
                    slascl_("G", &i_zero, &i_zero, &norm_a, &s_one, &m, &n, work, &ldwork, &INFO);
                }
                if(*trans == 'T')
                {
                    /*Copy x into work*/
                    fla_lapack_slacpy("All", &m, &nrhs, x, &ldb, &((float *)work)[n * ldwork],
                                      &ldwork);
                    norm_x = fla_lapack_slange("M", &m, &nrhs, &((float *)work)[n * ldwork],
                                               &ldwork, &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        slascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &m, &nrhs,
                                &((float *)work)[n * ldwork], &ldwork, &INFO);
                    }
                    temp = n + nrhs;
                    /*QR factorization of [A x]*/
                    sgeqr2_(&m, &temp, work, &ldwork, &((float *)work)[ldwork * (n + nrhs)],
                            &((float *)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))], &INFO);
                    norm = 0;
                    /*Compute norm*/
                    for(integer j = n + 1; j < temp; j++)
                    {
                        for(integer i = n + 1; i < fla_min(m, j); i++)
                        {
                            temp1 = FLA_FABS(((float *)work)[i + (j - 1) * m]);
                            norm = fla_max(temp1, norm);
                        }
                    }
                }
                else if(*trans == 'N')
                {
                    /*Copy x into work*/
                    for(int i = 0; i < n; i++)
                    {
                        for(int j = 0; j < nrhs; j++)
                        {
                            ((float *)work)[m + j + (i * ldwork)] = ((float *)x)[i + j * ldb];
                        }
                    }
                    norm_x
                        = fla_lapack_slange("M", &nrhs, &n, &((float *)work)[m], &ldwork, &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        slascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &nrhs, &n,
                                &((float *)work)[m + 1], &ldwork, &INFO);
                    }
                    /*LQ factorization [A' x]*/
                    sgelq2_(&ldwork, &n, work, &ldwork, &((float *)work)[ldwork * n],
                            &((float *)work)[ldwork * (n + 1)], &INFO);
                    /*Compute norm*/
                    for(integer j = m + 1; j < n; j++)
                    {
                        for(integer i = j; i < ldwork; i++)
                        {
                            temp1 = FLA_FABS(((float *)work)[i + (j * ldwork)]);
                            norm = fla_max(temp1, norm);
                        }
                    }
                }
                *(float *)resid = norm / ((double)fla_max(m, fla_max(n, nrhs)) * eps);
            }
            break;
        }
        case DOUBLE:
        {
            double eps, norm = 0, norm_a = 0, norm_x = 0, rwork;
            eps = fla_lapack_dlamch("E");

            if((*trans == 'T' && m > n) || (*trans == 'N' && m < n))
            {
                /*Copy A into work*/
                fla_lapack_dlacpy("All", &m, &n, A, &lda, work, &ldwork);
                norm_a = fla_lapack_dlange("M", &m, &n, work, &ldwork, &rwork);
                /*Scale work*/
                if(norm_a != 0.)
                {
                    dlascl_("G", &i_zero, &i_zero, &norm_a, &d_one, &m, &n, work, &ldwork, &INFO);
                }
                if(*trans == 'T')
                {
                    /*Copy x into work*/
                    fla_lapack_dlacpy("All", &m, &nrhs, x, &ldb, &((double *)work)[n * ldwork],
                                      &ldwork);
                    norm_x = fla_lapack_dlange("M", &m, &nrhs, &((double *)work)[n * ldwork],
                                               &ldwork, &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        dlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &m, &nrhs,
                                &((double *)work)[n * ldwork], &ldwork, &INFO);
                    }
                    temp = n + nrhs;
                    /*QR factorization of [A x]*/
                    dgeqr2_(&m, &temp, work, &ldwork, &((double *)work)[ldwork * (n + nrhs)],
                            &((double *)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))], &INFO);
                    norm = 0;
                    /*Compute norm*/
                    for(integer j = n + 1; j < temp; j++)
                    {
                        for(integer i = n + 1; i < fla_min(m, j); i++)
                        {
                            temp1 = FLA_FABS(((double *)work)[i + (j - 1) * m]);
                            norm = fla_max(temp1, norm);
                        }
                    }
                }
                else if(*trans == 'N')
                {
                    /*Copy x into work*/
                    for(int i = 0; i < n; i++)
                    {
                        for(int j = 0; j < nrhs; j++)
                        {
                            ((double *)work)[m + j + (i * ldwork)] = ((double *)x)[i + j * ldb];
                        }
                    }
                    norm_x
                        = fla_lapack_dlange("M", &nrhs, &n, &((double *)work)[m], &ldwork, &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        dlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &nrhs, &n,
                                &((double *)work)[m + 1], &ldwork, &INFO);
                    }
                    /*LQ factorization [A' x]*/
                    dgelq2_(&ldwork, &n, work, &ldwork, &((double *)work)[ldwork * n],
                            &((double *)work)[ldwork * (n + 1)], &INFO);
                    /*Compute norm*/
                    for(integer j = m + 1; j < n; j++)
                    {
                        for(integer i = j; i < ldwork; i++)
                        {
                            temp1 = FLA_FABS(((double *)work)[i + (j * ldwork)]);
                            norm = fla_max(temp1, norm);
                        }
                    }
                }
                *(double *)resid = norm / ((double)fla_max(m, fla_max(n, nrhs)) * eps);
            }
            break;
        }
        case COMPLEX:
        {
            float eps, norm = 0, norm_a = 0, norm_x = 0, rwork;
            eps = fla_lapack_slamch("E");
            if((*trans == 'C' && m > n) || (*trans == 'N' && m < n))
            {
                /*Copy A into work*/
                fla_lapack_clacpy("All", &m, &n, A, &lda, work, &ldwork);
                norm_a = fla_lapack_clange("M", &m, &n, work, &ldwork, &rwork);
                /*Scale work*/
                if(norm_a != 0.)
                {
                    clascl_("G", &i_zero, &i_zero, &norm_a, &s_one, &m, &n, work, &ldwork, &INFO);
                }
                if(*trans == 'C')
                {
                    /*Copy x into work*/
                    fla_lapack_clacpy("All", &m, &nrhs, x, &ldb, &((scomplex *)work)[n * ldwork],
                                      &ldwork);
                    norm_x = fla_lapack_clange("M", &m, &nrhs, &((scomplex *)work)[n * ldwork],
                                               &ldwork, &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        clascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &m, &nrhs,
                                &((scomplex *)work)[n * ldwork], &ldwork, &INFO);
                    }
                    temp = n + nrhs;
                    /*QR factorization of [A x]*/
                    cgeqr2_(&m, &temp, work, &ldwork, &((scomplex *)work)[ldwork * (n + nrhs)],
                            &((scomplex *)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))],
                            &INFO);
                    norm = 0;
                    /*Compute norm*/
                    for(integer j = n + 1; j < temp; j++)
                    {
                        for(integer i = n + 1; i < fla_min(m, j); i++)
                        {
                            temp1 = FLA_FABS(((scomplex *)work)[i + (j - 1) * m].real);
                            norm = fla_max(temp1, norm);
                        }
                    }
                }
                else if(*trans == 'N')
                {
                    /*Copy x into work*/
                    for(integer i = 0; i < n; i++)
                    {
                        for(integer j = 0; j < nrhs; j++)
                        {
                            ((scomplex *)work)[m + j + (i * ldwork)].real
                                = ((scomplex *)x)[i + j * ldb].real;
                            ((scomplex *)work)[m + j + (i * ldwork)].imag
                                = -1 * ((scomplex *)x)[i + j * ldb].imag;
                        }
                    }
                    norm_x = fla_lapack_clange("M", &nrhs, &n, &((scomplex *)work)[m], &ldwork,
                                               &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        clascl_("G", &i_zero, &i_zero, &norm_x, &s_one, &nrhs, &n,
                                &((scomplex *)work)[m + 1], &ldwork, &INFO);
                    }
                    /*LQ factorization [A' x]*/
                    cgelq2_(&ldwork, &n, work, &ldwork, &((scomplex *)work)[ldwork * n],
                            &((scomplex *)work)[ldwork * (n + 1)], &INFO);
                    /*Compute norm*/
                    for(integer j = m + 1; j < n; j++)
                    {
                        for(integer i = j; i < ldwork; i++)
                        {
                            temp1 = FLA_FABS(((scomplex *)work)[i + (j * ldwork)].real);
                            norm = fla_max(norm, temp1);
                        }
                    }
                }
                *(float *)resid = norm / ((double)fla_max(m, fla_max(n, nrhs)) * eps);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm = 0, norm_a = 0, norm_x = 0, rwork;
            eps = fla_lapack_dlamch("E");
            if((*trans == 'C' && m > n) || (*trans == 'N' && m < n))
            {
                /*Copy A into work*/
                fla_lapack_zlacpy("All", &m, &n, A, &lda, work, &ldwork);
                norm_a = fla_lapack_zlange("M", &m, &n, work, &ldwork, &rwork);
                /*Scale work*/
                if(norm_a != 0)
                {
                    zlascl_("G", &i_zero, &i_zero, &norm_a, &d_one, &m, &n, work, &ldwork, &INFO);
                }
                if(*trans == 'C')
                {
                    /*Copy x into work*/
                    fla_lapack_zlacpy("All", &m, &nrhs, x, &ldb, &((dcomplex *)work)[n * ldwork],
                                      &ldwork);
                    norm_x = fla_lapack_zlange("M", &m, &nrhs, &((dcomplex *)work)[n * ldwork],
                                               &ldwork, &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        zlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &m, &nrhs,
                                &((dcomplex *)work)[n * ldwork], &ldwork, &INFO);
                    }
                    temp = n + nrhs;
                    /*QR factorization of [A x]*/
                    zgeqr2_(&m, &temp, work, &ldwork, &((dcomplex *)work)[ldwork * (n + nrhs)],
                            &((dcomplex *)work)[ldwork * (n + nrhs) + fla_min(m, (n + nrhs))],
                            &INFO);
                    norm = 0;
                    /*Compute norm*/
                    for(integer j = n + 1; j < temp; j++)
                    {
                        for(integer i = n + 1; i < fla_min(m, j); i++)
                        {
                            temp1 = FLA_FABS(((dcomplex *)work)[i + (j - 1) * m].real);
                            norm = fla_max(temp1, norm);
                        }
                    }
                }
                else if(*trans == 'N')
                {
                    /*Copy x into work*/
                    for(integer i = 0; i < n; i++)
                    {
                        for(integer j = 0; j < nrhs; j++)
                        {
                            ((dcomplex *)work)[m + j + (i * ldwork)] = ((dcomplex *)x)[i + j * ldb];
                            ((dcomplex *)work)[m + j + (i * ldwork)].imag
                                = -1 * ((dcomplex *)x)[i + j * ldb].imag;
                        }
                    }
                    norm_x = fla_lapack_zlange("M", &nrhs, &n, &((dcomplex *)work)[m], &ldwork,
                                               &rwork);
                    /*Scale x*/
                    if(norm_x != 0)
                    {
                        zlascl_("G", &i_zero, &i_zero, &norm_x, &d_one, &nrhs, &n,
                                &((dcomplex *)work)[m + 1], &ldwork, &INFO);
                    }
                    /*LQ factorization [A' x]*/
                    zgelq2_(&ldwork, &n, work, &ldwork, &((dcomplex *)work)[ldwork * n],
                            &((dcomplex *)work)[ldwork * (n + 1)], &INFO);
                    /*Compute norm*/
                    for(integer j = m + 1; j < n; j++)
                    {
                        for(integer i = j; i < ldwork; i++)
                        {
                            temp1 = FLA_FABS(((dcomplex *)work)[i + (j * ldwork)].real);
                            norm = fla_max(norm, temp1);
                        }
                    }
                }
                *(double *)resid = norm / ((double)fla_max(m, fla_max(n, nrhs)) * eps);
            }
            break;
        }
    }
    free_vector(work);
}

/* To calculate the resudial sum of squares of solution for solution x of Ax = b and m < n
 */
void residual_sum_of_squares(int datatype, integer m, integer n, integer nrhs, void *x, integer ldx,
                             double *resid)
{
    integer temp = m - n;
    *resid = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            for(integer i = 0; i < nrhs; i++)
            {
                *resid = fla_max(snrm2_(&temp, &((float *)x)[(i * ldx) + n], &i_one), *resid);
            }
            break;
        }
        case DOUBLE:
        {
            for(integer i = 0; i < nrhs; i++)
            {
                *resid = fla_max(dnrm2_(&temp, &((double *)x)[(i * ldx) + n], &i_one), *resid);
            }
            break;
        }
        case COMPLEX:
        {
            for(integer i = 0; i < nrhs; i++)
            {
                *resid = fla_max(scnrm2_(&temp, &((scomplex *)x)[(i * ldx) + n], &i_one), *resid);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(integer i = 0; i < nrhs; i++)
            {
                *resid = fla_max(dznrm2_(&temp, &((dcomplex *)x)[(i * ldx) + n], &i_one), *resid);
            }
            break;
        }
    }
}

/* swap row or column in a matrix
 * This function can be used to swap -> a row of the matrix with another row
                                     -> a column of matrix with another column
                                     -> a row with a column and vice-versa
 * m - Elements in vector to be copied
 * A - Matrix
 * incx - Increment for first row/col
 * incy - Increment for second row/col
 * (srow, scol) - Start location of the first vector in a matrix
 * (drow, dcol) - Start location of the second vector in a matrix */
void swap_row_col(integer datatype, integer *m, void *A, integer lda, integer *incx, integer *incy,
                  integer srow, integer scol, integer drow, integer dcol)
{
    void *x = NULL, *y = NULL;
    switch(datatype)
    {
        case FLOAT:
        {
            x = ((float *)A + (scol * lda + srow));
            y = ((float *)A + (dcol * lda + drow));
            sswap_(m, x, incx, y, incy);
            break;
        }
        case DOUBLE:
        {
            x = ((double *)A + (scol * lda + srow));
            y = ((double *)A + (dcol * lda + drow));
            dswap_(m, x, incx, y, incy);
            break;
        }
        case COMPLEX:
        {
            x = ((scomplex *)A + (scol * lda + srow));
            y = ((scomplex *)A + (dcol * lda + drow));
            cswap_(m, x, incx, y, incy);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            x = ((dcomplex *)A + (scol * lda + srow));
            y = ((dcomplex *)A + (dcol * lda + drow));
            zswap_(m, x, incx, y, incy);
            break;
        }
    }
}

/* GEMM implementation for  C := alpha*op( A )*op( B ) + beta*C
 * Where alpha = 1, beta = 0
 */
void fla_invoke_gemm(integer datatype, char *transA, char *transB, integer *m, integer *n,
                     integer *k, void *A, integer *lda, void *B, integer *ldb, void *C,
                     integer *ldc)
{
    switch(datatype)
    {
        case FLOAT:
        {
            sgemm_(transA, transB, m, n, k, &s_one, A, lda, B, ldb, &s_zero, C, ldc);
            break;
        }
        case DOUBLE:
        {
            dgemm_(transA, transB, m, n, k, &d_one, A, lda, B, ldb, &d_zero, C, ldc);
            break;
        }
        case COMPLEX:
        {
            cgemm_(transA, transB, m, n, k, &c_one, A, lda, B, ldb, &c_zero, C, ldc);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zgemm_(transA, transB, m, n, k, &z_one, A, lda, B, ldb, &z_zero, C, ldc);
            break;
        }
    }
}

/* Generate a symmetric or hermitian matrix from existing matrix A
 * If type = "C" hermitian matrix formed.
 * If type = "S" symmetric matrix is formed.
 */
void form_symmetric_matrix(integer datatype, integer n, void *A, integer lda, char *type)
{
    integer i, j, conj = 1;
    if(*type == 'C')
    {
        conj = -1;
    }
    switch(datatype)
    {
        case FLOAT:
        {
            for(i = 0; i < n; i++)
            {
                for(j = i; j < n; j++)
                {
                    ((float *)A)[i * lda + j] = ((float *)A)[j * lda + i];
                }
            }
            break;
        }
        case DOUBLE:
        {
            for(i = 0; i < n; i++)
            {
                for(j = i; j < n; j++)
                {
                    ((double *)A)[i * lda + j] = ((double *)A)[j * lda + i];
                }
            }
            break;
        }
        case COMPLEX:
        {
            for(i = 0; i < n; i++)
            {
                for(j = i; j < n; j++)
                {
                    if(i == j)
                    {
                        continue;
                    }
                    ((scomplex *)A)[j * lda + i].real = ((scomplex *)A)[i * lda + j].real;
                    ((scomplex *)A)[j * lda + i].imag = conj * ((scomplex *)A)[i * lda + j].imag;
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for(i = 0; i < n; i++)
            {
                for(j = i; j < n; j++)
                {
                    if(i == j)
                    {
                        continue;
                    }
                    ((dcomplex *)A)[j * lda + i].real = ((dcomplex *)A)[i * lda + j].real;
                    ((dcomplex *)A)[j * lda + i].imag = conj * ((dcomplex *)A)[i * lda + j].imag;
                }
            }
            break;
        }
    }
}
/* Scaling the matrix by x scalar */
void scal_matrix(integer datatype, void *x, void *A, integer m, integer n, integer lda, integer inc)
{
    integer j;
    switch(datatype)
    {
        case FLOAT:
        {
            float *column_of_matrix;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by scalar x */
                column_of_matrix = (float *)A + j * lda;
                sscal_(&m, x, column_of_matrix, &inc);
            }
            break;
        }

        case DOUBLE:
        {
            double *column_of_matrix;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by in the x */
                column_of_matrix = (double *)A + j * lda;
                dscal_(&m, x, column_of_matrix, &inc);
            }
            break;
        }

        case COMPLEX:
        {
            scomplex *column_of_matrix;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by x */
                column_of_matrix = (scomplex *)A + j * lda;
                csscal_(&m, x, column_of_matrix, &inc);
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            dcomplex *column_of_matrix;
            for(j = 0; j < n; j++)
            {
                /* scale each column of the matrix by x */
                column_of_matrix = (dcomplex *)A + j * lda;
                zdscal_(&m, x, column_of_matrix, &inc);
            }
            break;
        }
    }
    return;
}
/* Get the maximum value from the matrix */
void get_max_from_matrix(integer datatype, void *A, void *max_val, integer m, integer n,
                         integer lda)
{
    void *work = NULL;
    switch(datatype)
    {
        case FLOAT:
        {
            *(float *)max_val = fla_lapack_slange("M", &m, &n, A, &lda, work);
            break;
        }

        case DOUBLE:
        {
            *(double *)max_val = fla_lapack_dlange("M", &m, &n, A, &lda, work);
            break;
        }
        case COMPLEX:
        {
            *(float *)max_val = fla_lapack_clange("M", &m, &n, A, &lda, work);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            *(double *)max_val = fla_lapack_zlange("M", &m, &n, A, &lda, work);
            break;
        }
    }
}

/* Get the minimum value from the matrix */
void get_min_from_matrix(integer datatype, void *A, void *min_val, integer m, integer n,
                         integer lda)
{
    integer i, j;

    switch(datatype)
    {
        case FLOAT:
        {
            float min_local = FLT_MAX;
            for(i = 0; i < n; i++)
            {
                for(j = 0; j < m; j++)
                {
                    float val = FLA_FABS(((float *)A)[i * lda + j]);
                    if(val != s_zero && min_local > val)
                    {
                        min_local = val;
                    }
                }
            }

            *(float *)min_val = min_local;
            break;
        }

        case DOUBLE:
        {
            double min_local = DBL_MAX;
            for(i = 0; i < n; i++)
            {
                for(j = 0; j < m; j++)
                {
                    double val = FLA_FABS(((double *)A)[i * lda + j]);
                    if(val != d_zero && min_local > val)
                    {
                        min_local = val;
                    }
                }
            }
            *(double *)min_val = min_local;
            break;
        }

        case COMPLEX:
        {
            scomplex *ptr = A;
            scomplex min_local;
            min_local.real = FLT_MAX;
            min_local.imag = FLT_MAX;

            for(i = 0; i < n; i++)
            {
                for(j = 0; j < m; j++)
                {
                    float real, imag;
                    real = FLA_FABS(ptr[i * lda + j].real);
                    imag = FLA_FABS(ptr[i * lda + j].imag);
                    /* Compare real part */
                    if(real != s_zero && real < min_local.real)
                    {
                        min_local.real = real;
                    }
                    /* Compare imaginary part */
                    if(imag != s_zero && imag < min_local.imag)
                    {
                        min_local.imag = imag;
                    }
                }
            }
            *(float *)min_val = fla_min(min_local.real, min_local.imag);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            dcomplex *ptr = A;
            dcomplex min_local;
            min_local.real = DBL_MAX;
            min_local.imag = DBL_MAX;

            for(i = 0; i < n; i++)
            {
                for(j = 0; j < m; j++)
                {
                    double real, imag;
                    real = FLA_FABS(ptr[i * lda + j].real);
                    imag = FLA_FABS(ptr[i * lda + j].imag);
                    /* Compare real part */
                    if(real != d_zero && real < min_local.real)
                    {
                        min_local.real = real;
                    }
                    /* Compare imaginary part */
                    if(imag != d_zero && imag < min_local.imag)
                    {
                        min_local.imag = imag;
                    }
                }
            }
            *(double *)min_val = fla_min(min_local.real, min_local.imag);
            break;
        }
    }
}

/* Sort the input vector in the given order
   order = A - Ascending order
         = D - Descending order */
void sort_vector(integer datatype, char *order, integer vect_len, void *w, integer incw)
{
    if(!w)
        return;

    integer i, j;
    switch(datatype)
    {
        case FLOAT:
        {
            float temp;
            float *w_ptr = (float *)w;

            for(i = 0; i < vect_len; i++)
            {
                for(j = i + 1; j < vect_len; j++)
                {
                    if(*order == 'A')
                    {
                        if(*(w_ptr + i * incw) > *(w_ptr + j * incw))
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                    }
                    else if(*order == 'D')
                    {
                        if(*(w_ptr + i * incw) < *(w_ptr + j * incw))
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {

            double temp;
            double *w_ptr = (double *)w;

            for(i = 0; i < vect_len; i++)
            {
                for(j = i + 1; j < vect_len; j++)
                {
                    if(*order == 'A')
                    {
                        if(*(w_ptr + i * incw) > *(w_ptr + j * incw))
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                    }
                    else if(*order == 'D')
                    {
                        if(*(w_ptr + i * incw) < *(w_ptr + j * incw))
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            scomplex temp;
            scomplex *w_ptr = (scomplex *)w;

            for(i = 0; i < vect_len; i++)
            {
                for(j = i + 1; j < vect_len; j++)
                {
                    if(*order == 'A')
                    {
                        if((w_ptr + i * incw)->real > (w_ptr + j * incw)->real)
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                        else if((w_ptr + i * incw)->real == (w_ptr + j * incw)->real)
                        {
                            if((w_ptr + i * incw)->imag > (w_ptr + j * incw)->imag)
                            {
                                temp = *(w_ptr + i * incw);
                                *(w_ptr + i * incw) = *(w_ptr + j * incw);
                                *(w_ptr + j * incw) = temp;
                            }
                        }
                    }
                    else if(*order == 'D')
                    {
                        if((w_ptr + i * incw)->real < (w_ptr + j * incw)->real)
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                        else if((w_ptr + i * incw)->real == (w_ptr + j * incw)->real)
                        {
                            if((w_ptr + i * incw)->imag < (w_ptr + j * incw)->imag)
                            {
                                temp = *(w_ptr + i * incw);
                                *(w_ptr + i * incw) = *(w_ptr + j * incw);
                                *(w_ptr + j * incw) = temp;
                            }
                        }
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex temp;
            dcomplex *w_ptr = (dcomplex *)w;

            for(i = 0; i < vect_len; i++)
            {
                for(j = i + 1; j < vect_len; j++)
                {
                    if(*order == 'A')
                    {
                        if((w_ptr + i * incw)->real > (w_ptr + j * incw)->real)
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                        else if((w_ptr + i * incw)->real == (w_ptr + j * incw)->real)
                        {
                            if((w_ptr + i * incw)->imag > (w_ptr + j * incw)->imag)
                            {
                                temp = *(w_ptr + i * incw);
                                *(w_ptr + i * incw) = *(w_ptr + j * incw);
                                *(w_ptr + j * incw) = temp;
                            }
                        }
                    }
                    else if(*order == 'D')
                    {
                        if((w_ptr + i * incw)->real < (w_ptr + j * incw)->real)
                        {
                            temp = *(w_ptr + i * incw);
                            *(w_ptr + i * incw) = *(w_ptr + j * incw);
                            *(w_ptr + j * incw) = temp;
                        }
                        else if((w_ptr + i * incw)->real == (w_ptr + j * incw)->real)
                        {
                            if((w_ptr + i * incw)->imag < (w_ptr + j * incw)->imag)
                            {
                                temp = *(w_ptr + i * incw);
                                *(w_ptr + i * incw) = *(w_ptr + j * incw);
                                *(w_ptr + j * incw) = temp;
                            }
                        }
                    }
                }
            }
            break;
        }
    }
}

/* Generate a block diagonal matrix with complex conjugate eigen value pairs as
   2 * 2 blocks along the diagonal. This is used for generating asymmetric matrix */
void create_realtype_block_diagonal_matrix(integer datatype, void *A, integer n, integer lda)
{
    integer i;
    if(lda < n)
        return;

    reset_matrix(datatype, n, n, A, lda);

    if(datatype == FLOAT)
    {
        for(i = 0; i < n; i += 2)
        {
            ((float *)A)[i * lda + i] = SRAND();
            if(i < n-1)
            {
                ((float *)A)[(i + 1) * lda + (i + 1)] = ((float *)A)[i * lda + i];
                ((float *)A)[i * lda + i + 1] = SRAND();
                ((float *)A)[(i + 1) * lda + i] = -((float *)A)[i * lda + i + 1];
            }
        }
    }
    else if(datatype == DOUBLE)
    {
        for(i = 0; i < n; i += 2)
        {
            ((double *)A)[i * lda + i] = DRAND();
            if(i < n-1)
            {
                ((double *)A)[(i + 1) * lda + (i + 1)] = ((double *)A)[i * lda + i];
                ((double *)A)[i * lda + i + 1] = DRAND();
                ((double *)A)[(i + 1) * lda + i] = -((double *)A)[i * lda + i + 1];
            }
        }
    }
}

/*
 *   Create input matrix A by randomly generating eigen values(EVs)
 *               A  = (Q * lambda * Q')
 *   where  lambda  is a super diagonal matrix with diagonal, sub diagonal elements being
 *                  the eigen values of input matrix A
 *               Q  is an orthogonal matrix with corresponding
 *                  eigen vectors as its rows.
 */
void generate_asym_matrix_from_EVs(integer datatype, integer n, void *A, integer lda, void *L,
                                   char imatrix, void *scal)
{
    void *X = NULL, *Q = NULL;
    integer info = 0;

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        void *L1 = NULL;
        create_vector(datatype, &L1, n);
        rand_vector(datatype, n, L1, 1, d_zero, d_zero, 'R');
        diagonalize_vector(datatype, L1, L, n, n, n);

        free_vector(L1);
    }
    else
    {
        create_realtype_block_diagonal_matrix(datatype, L, n, n);
    }

    create_matrix(datatype, &X, n, n);
    rand_matrix(datatype, X, n, n, n);
    create_matrix(datatype, &Q, n, n);

    /* Generate random orthogonal matrix(Q) of size n x n */
    get_orthogonal_matrix_from_QR(datatype, n, X, n, Q, n, &info);

    /* Generate input matrix A using L(Eigen values)
       and Q(Eigen vectors) obtained above using reverse
       Eigen decompostion */
    generate_asym_matrix_from_ED(datatype, n, A, lda, Q, L);

    if(imatrix == 'O' || imatrix == 'U')
        init_matrix_overflow_underflow_asym(datatype, n, n, A, lda, imatrix, scal);

    /* Free up the buffers */
    free_matrix(X);
    free_matrix(Q);
    return;
}

/*
 * Generate square matrix of size n x n using Eigen decomposition(ED)
 *                     A  = (Q * lambda * Q')
 * where Q is an n x n orthogonal matrix and
 *       lambda is a block diagonal matrix or triangular matrix of size n * n
 */
void generate_asym_matrix_from_ED(integer datatype, integer n, void *A, integer lda, void *Q,
                                  void *lambda)
{
    void *Qlambda = NULL;
    if(lda < n)
        return;

    create_matrix(datatype, &Qlambda, n, n);

    switch(datatype)
    {
        case FLOAT:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            sgemm_("N", "N", &n, &n, &n, &s_one, Q, &n, lambda, &n, &s_zero, Qlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, Qlambda, &n, Q, &n, &s_zero, A, &lda);
            break;
        }
        case DOUBLE:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            dgemm_("N", "N", &n, &n, &n, &d_one, Q, &n, lambda, &n, &d_zero, Qlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, Qlambda, &n, Q, &n, &d_zero, A, &lda);
            break;
        }
        case COMPLEX:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            cgemm_("N", "N", &n, &n, &n, &c_one, Q, &n, lambda, &n, &c_zero, Qlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, Qlambda, &n, Q, &n, &c_zero, A, &lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Generate matrix A using eigen decomposition (Q * lambda) * Q' */
            zgemm_("N", "N", &n, &n, &n, &z_one, Q, &n, lambda, &n, &z_zero, Qlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, Qlambda, &n, Q, &n, &z_zero, A, &lda);
            break;
        }
    }

    free_matrix(Qlambda);
}

/* Compare two vectors starting from offset_A in A vector with B vector
   (starting from offset 0 in B) */
integer compare_vector(integer datatype, integer vect_len, void *A, integer inca, integer offset_A,
                       void *B, integer incb)
{
    integer i;
    switch(datatype)
    {
        case FLOAT:
        {
            float *a = (float *)A;
            float *b = (float *)B;
            for(i = 0; i < vect_len; i++)
            {
                if(f2c_abs(a[(i * inca) + (offset_A - 1)] - b[i * incb]) > MAX_FLT_DIFF)
                {
                    return 1;
                }
            }
            break;
        }
        case DOUBLE:
        {
            double *a = (double *)A;
            double *b = (double *)B;
            for(i = 0; i < vect_len; i++)
            {
                if(f2c_abs(a[(i * inca) + (offset_A - 1)] - b[i * incb]) > MAX_DBL_DIFF)
                {
                    return 1;
                }
            }
            break;
        }
        case COMPLEX:
        {
            scomplex *a = (scomplex *)A;
            scomplex *b = (scomplex *)B;
            for(i = 0; i < vect_len; i++)
            {
                if(f2c_abs(a[(i * inca) + (offset_A - 1)].real - b[i * incb].real) > MAX_FLT_DIFF)
                {
                    return 1;
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex *a = (dcomplex *)A;
            dcomplex *b = (dcomplex *)B;
            for(i = 0; i < vect_len; i++)
            {
                if(f2c_abs(a[(i * inca) + (offset_A - 1)].real - b[i * incb].real) > MAX_DBL_DIFF)
                {
                    return 1;
                }
            }
            break;
        }
    }
    return 0;
}

/* Create diagonal matrix by copying elements from a vector to matrix */
void diagonalize_vector(integer datatype, void *s, void *sigma, integer m, integer n, integer LDA)
{
    integer incr, min_m_n;

    incr = m + 1;
    min_m_n = fla_min(m, n);

    reset_matrix(datatype, m, n, sigma, m);

    switch(datatype)
    {
        case FLOAT:
        {
            scopy_(&min_m_n, s, &i_one, sigma, &incr);
            break;
        }
        case DOUBLE:
        {
            dcopy_(&min_m_n, s, &i_one, sigma, &incr);
            break;
        }
        case COMPLEX:
        {
            ccopy_(&min_m_n, s, &i_one, sigma, &incr);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zcopy_(&min_m_n, s, &i_one, sigma, &incr);
            break;
        }
    }
}

/* Find negative value of each element and store in next location
   Used to store imaginary parts of complex conjuate pair of eigen values
   in asymmetric matrix eigen decomposition APIs
   Ex: input vector {a, 0, -b, 0 ...}
       output vector {a, -a, -b, b, ...} */
void add_negative_values(integer datatype, void *vect, integer n)
{
    if(!vect)
        return;

    integer i;
    if(datatype == FLOAT)
    {
        float *w = (float *)vect;
        for(i = 0; i < n - 1; i++)
        {
            if(i % 2 == 0)
            {
                w[i + 1] = -w[i];
            }
        }
    }
    else
    {
        double *w = (double *)vect;
        for(i = 0; i < n - 1; i++)
        {
            if(i % 2 == 0)
            {
                w[i + 1] = -w[i];
            }
        }
    }
}
