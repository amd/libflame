/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_steqr_experiment(test_params_t *params, integer datatype, integer p_cur, integer  q_cur, integer pci,
integer n_repeats, double* perf, double* t, double* residual);
void prepare_steqr_run(char* compz, integer n, void* Z, void* D, void* E, integer datatype, integer n_repeats, double* time_min_);
void invoke_steqr(integer datatype, char* compz, integer* n, void* z, integer* ldz, void* d, void* e, void* work, integer* info);

void fla_test_steqr(test_params_t *params)
{
    char* op_str = "Eigen Decomposition of symmetrix tridiagonal matrix";
    char* front_str = "STEQR";

    fla_test_output_info("--- %s ---\n", op_str);
    fla_test_output_info("\n");
    fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_steqr_experiment);
}

void fla_test_steqr_experiment(test_params_t *params,
                               integer  datatype,
                               integer  p_cur,
                               integer  q_cur,
                               integer pci,
                               integer n_repeats,
                               double* perf,
                               double *time_min,
                               double* residual)
{
    integer n, ldz, lda, info = 0;
    char compz, uplo;
    void *Z = NULL, *Z_test = NULL, *A = NULL, *Q = NULL;
    void *D = NULL, *D_test = NULL, *E = NULL, *E_test = NULL;

    /* Get input matrix dimensions.*/
    compz = params->eig_sym_paramslist[pci].compz;
    *residual = params->eig_sym_paramslist[pci].threshold_value;
    uplo = params->eig_sym_paramslist[pci].uplo;

    n = p_cur;
    ldz = max(1,n);
    lda = max(1,n);

    /* Create input matrix parameters */
    create_matrix(datatype, &Z, n, n);
    create_matrix(datatype, &A, n, n);
    create_matrix(datatype, &Q, n, n);

    reset_matrix(datatype, n, n, Z, ldz);
    reset_matrix(datatype, n, n, A, lda);
    reset_matrix(datatype, n, n, Q, ldz);

    create_vector(get_realtype(datatype), &D, n);
    create_vector(get_realtype(datatype), &E, n-1);

    /* input matrix Z with random symmetric numbers and D,E matrix with diagonal and subdiagonal values */
    if(datatype == FLOAT || datatype == DOUBLE)
        rand_sym_matrix(datatype, A, n, n, lda);
    else
        rand_hermitian_matrix(datatype, n, &A, lda);

    copy_matrix(datatype, "full", n, n, A, lda, Q, ldz);
    /* Make a copy of input matrix Z. This is required to validate the API functionality.*/
    create_matrix(datatype, &Z_test, n, n);
    reset_matrix(datatype, n, n, Z_test, ldz);

    invoke_sytrd(datatype, &uplo, compz, n, Q, lda, D, E, info);
    /*form tridiagonal matrix Z by copying from matrix*/
    copy_sym_tridiag_matrix(datatype, D, E, n, n, Z, ldz);

    if(compz == 'I')
        set_identity_matrix(datatype, n, n, Z_test, ldz);
    else
    {
        copy_matrix(datatype, "full", n, n, Q, lda, Z_test, ldz);
        copy_matrix(datatype, "full", n, n, A, lda, Z, ldz);
    }
    create_vector(get_realtype(datatype), &D_test, n);
    create_vector(get_realtype(datatype), &E_test, n-1);
    copy_vector(get_realtype(datatype), n, D, 1, D_test, 1);
    copy_vector(get_realtype(datatype), n-1, E, 1, E_test, 1);

    prepare_steqr_run(&compz, n, Z_test, D_test, E_test, datatype, n_repeats, time_min);

    /* performance computation
       24 n^2 flops for eigen vectors of Z, compz = 'N'
       7 n^3 flops for eigen vectors of Z, compz = 'V' or 'I'
       14 n^3 flops for eigen vectors of Z for complex, compz = 'V' or 'I' */
 
    if( compz == 'I' || compz == 'V')
        *perf = (double)(7.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    else if( compz == 'N')
        *perf = (double)(24.0 * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf = (double)(14.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;

    /* output validation */
    validate_syevd(&compz, n, Z, Z_test, D_test, datatype, residual);

    /* Free up the buffers */
    free_matrix(Z);
    free_vector(D);
    free_vector(E);
    free_matrix(A);
    free_matrix(Q);
    free_matrix(Z_test);
    free_vector(D_test);
    free_vector(E_test);
}

void prepare_steqr_run(char *compz,
                       integer n,
                       void *Z,
                       void *D,
                       void *E,
                       integer datatype,
                       integer n_repeats,
                       double* time_min_)
{
    integer ldz;
    void *Z_save = NULL, *D_save = NULL, *E_save = NULL, *work = NULL;
    integer i, info = 0;
    double time_min = 1e9, exe_time;

    ldz = max(1,n);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &Z_save, n, n);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);

    create_vector(get_realtype(datatype), &D_save, n);
    create_vector(get_realtype(datatype), &E_save, n-1);
    copy_vector(get_realtype(datatype), n, D, 1, D_save, 1);
    copy_vector(get_realtype(datatype), n-1, E, 1, E_save, 1);

    for (i = 0; i < n_repeats; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        copy_vector(get_realtype(datatype), n, D_save, 1, D, 1);
        copy_vector(get_realtype(datatype), n-1, E_save, 1, E, 1);

        create_vector(get_realtype(datatype), &work, 2 * (n - 1));
        exe_time = fla_test_clock();

        /* call to API */
        invoke_steqr(datatype, compz, &n, Z, &ldz, D, E, work, &info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }

    *time_min_ = time_min;

    free(Z_save);
    free(D_save);
    free(E_save);
}

void invoke_steqr(integer datatype, char* compz, integer* n, void* z, integer* ldz, void* d, void* e, void* work, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            ssteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case DOUBLE:
        {
            dsteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case COMPLEX:
        {
            csteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zsteqr_(compz, n, d, e, z, ldz, work, info);
            break;
        }
    }
}