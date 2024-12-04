/*
    Copyright (C) 2023 - 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*
    Computes the LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.
*/

#include "FLAME.h"

#include <stdio.h>
#include <stdlib.h>

/* Generate Random Value */
#define DRAND() ((double)rand() / ((double)RAND_MAX / 2.0F)) - 1.0F;

/* Local Function Declaration */
void rand_matrix(void *A, integer M, integer N, integer LDA);

int main(int argc, char **argv)
{
    /* Initialize input matrix sizes */
    integer M = 10, N = 10, LDA = 10;
    double *A;
    integer *ipiv;
    integer info;

    /* Allocation of memory to matrix*/
    A = (double *)malloc(LDA * N * sizeof(double));
    ipiv = (integer *)malloc(N * sizeof(double));

    /* Initialize matrix with random values */
    rand_matrix(A, M, N, LDA);

    printf("Started execution of DGETRF API \n");

    /* Call to the DGETRF API */
    dgetrf_(&M, &N, A, &LDA, ipiv, &info);

    if(info == 0)
        printf("DGETRF API execution sucessfully completed\n");
    else
        printf("DGETRF Execution Failed\n");
    return 0;
}

/* Initialize the matrix with random values */
void rand_matrix(void *A, integer M, integer N, integer LDA)
{
    integer i, j;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < M; j++)
        {
            ((double *)A)[i * LDA + j] = DRAND();
        }
    }
}
