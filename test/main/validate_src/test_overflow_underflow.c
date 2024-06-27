/*
    Copyright (C) 2024, Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "test_overflow_underflow.h"

/*
 * Step 1: Create matrix with properties specific to API
 * Step 2: Find ratio between the  max, min value for datatype
 *        and matrix max min value.
 * Step 3: Scale the matrix by x scalar
 * Explanation w.r.t overflow amd underflow
 * In overflow cases :
 * Calculate the scaling value as ratio between
 * maximum value for datatype and maximum value within the matrix
 *
 * FLT_MAX = 3.402823E+38
 * DBL_MAX = 1.797693E+308
 *
 *       *************************************
 *       * RATIO = DATATYPE_MAX / MATRIX_MAX *
 *       *************************************
 * Scale up  the matrix with ratio as scaling factor
 *       scale = ratio
 *       *************************************
 *       *       A' =  ( scale * A )         *
 *       *************************************
 *
 * In underflow cases :
 * Calculate the scaling value as ratio between
 * minimum value for datatype and minimum value within the matrix
 *
 * FLT_MIN = 1.175494E-38
 * DBL_MIN = 2.225074E-308
 *
 *       *************************************
 *       * RATIO = DATATYPE_MIN / MATRIX_MIN *
 *       *************************************
 * Scale down the matrix with ratio as scaling factor
 *       scale = ratio
 *       *************************************
 *       *       A' =  ( scale * A )         *
 *       *************************************
 */

/* Initializing svd matrix with values around overflow underflow (GESVD) */
void init_matrix_overflow_underflow_svd(integer datatype, integer m, integer n, void *A,
                                        integer lda, char imatrix, void *scal)
{
    /* Find factor to scale up/down the matrix */
    calculate_svd_scale_value(datatype, m, n, A, lda, imatrix, scal);
    /* Scaling the matrix A by x scalar */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);
}

/* Calculating the scaling value with respect to max and min for SVD */
void calculate_svd_scale_value(integer datatype, integer m, integer n, void *A, integer lda,
                               char imatrix, void *scal)
{
    float flt_ratio;
    double dbl_ratio;

    if(imatrix == 'O')
    {
        /*Intialize the ratios with maximum of datatype value*/
        flt_ratio = FLT_MAX;
        dbl_ratio = DBL_MAX;
        void *max;
        create_vector(get_realtype(datatype), &max, 1);
        get_max_from_matrix(datatype, A, max, m, n, lda);
        /* The ratio is modified w.r.t dimension, to avoid inf in output */
        if(fla_max(m, n) <= 50)
        {
            flt_ratio = flt_ratio / 7.00;
            dbl_ratio = dbl_ratio / 10.00;
        }
        else if(fla_max(m, n) <= 500)
        {
            flt_ratio = flt_ratio / 11.00;
            dbl_ratio = dbl_ratio / 10.00;
        }
        else if(fla_max(m, n) <= 1000)
        {
            flt_ratio = flt_ratio / 15.00;
            dbl_ratio = dbl_ratio / 17.00;
        }
        else
        {
            flt_ratio = flt_ratio / 50.00;
            dbl_ratio = dbl_ratio / 50.00;
        }
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, max);
        free_vector(max);
    }
    else if(imatrix == 'U')
    {
        void *min;
        create_vector(get_realtype(datatype), &min, 1);
        /* Get minimum value from matrix */
        get_min_from_matrix(datatype, A, min, m, n, lda);
        /* Intialize with minimum of datatype value */
        flt_ratio = FLT_MIN;
        dbl_ratio = DBL_MIN;
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, min);
        free_vector(min);
    }
}

/* Initializing asymmetrix matrix with values around overflow underflow (GEEV)*/
void init_matrix_overflow_underflow_asym(integer datatype, integer m, integer n, void *A,
                                         integer lda, char imatrix, void *scal)
{
    /* Find factor to scale up/down the matrix */
    calculate_asym_scale_value(datatype, m, n, A, lda, imatrix, scal);
    /* Scaling the matrix A by X scalar */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);
}

/* Calculating the scaling value with respect to max and min for ASYM */
void calculate_asym_scale_value(integer datatype, integer m, integer n, void *A, integer lda,
                                char imatrix, void *scal)
{
    float flt_ratio;
    double dbl_ratio;
    if(n == 1)
    {
        assign_value(datatype, A, d_one, d_one);
    }
    if(imatrix == 'O')
    {
        /*Intialize the ratios with maximum of datatype value*/
        flt_ratio = FLT_MAX;
        dbl_ratio = DBL_MAX;
        void *max;
        create_vector(get_realtype(datatype), &max, 1);
        get_max_from_matrix(datatype, A, max, m, n, lda);
        /* The ratio is modified w.r.t dimension, to avoid inf in output */
        if(n <= 50)
        {
            flt_ratio = FLT_MAX / 4.00;
            dbl_ratio = DBL_MAX / 4.00;
        }
        else if(n <= 100)
        {
            flt_ratio = flt_ratio / 6.00;
            dbl_ratio = dbl_ratio / 6.00;
        }
        else if(n <= 500)
        {
            flt_ratio = flt_ratio / 15.00;
            dbl_ratio = dbl_ratio / 15.00;
        }
        else if(n <= 1000)
        {
            flt_ratio = flt_ratio / 20.00;
            dbl_ratio = dbl_ratio / 20.00;
        }
        else
        {
            flt_ratio = flt_ratio / 1000.00;
            dbl_ratio = dbl_ratio / 1000.00;
        }
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, max);
        free_vector(max);
    }
    else if(imatrix == 'U')
    {
        void *min;
        create_vector(get_realtype(datatype), &min, 1);
        /* Get minimum value from matrix */
        get_min_from_matrix(datatype, A, min, m, n, lda);
        /* Intialize with minimum of datatype value */
        flt_ratio = FLT_MIN;
        dbl_ratio = DBL_MIN;
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, min);
        free_vector(min);
    }
}

/* Finding ratio between two real datatype values */
void compute_ratio(integer datatype, void *scal, float flt_quotient, double dbl_quotient,
                   void *denominator)
{
    if(datatype == FLOAT)
    {
        *(float *)scal = flt_quotient / *(float *)denominator;
    }
    else if(datatype == DOUBLE)
    {
        *(double *)scal = dbl_quotient / *(double *)denominator;
    }
}

/* calculate the scaling value based on datatype max / min, matrix max val/ min val and a tuning
 * value */
void calculate_scale_value(integer datatype, void *scal, void *max_min, double tuning_val,
                           char imatrix_char)
{
    float feasible_float = FLT_MAX;
    double feasible_double = DBL_MAX;
    if(imatrix_char == 'U')
    {
        feasible_float = FLT_MIN;
        feasible_double = DBL_MIN;
    }

    if(get_realtype(datatype) == FLOAT)
    {
        *(float *)scal = (feasible_float / tuning_val) / *(float *)max_min;
    }
    else if(get_realtype(datatype) == DOUBLE)
    {
        *(double *)scal = (feasible_double / tuning_val) / *(double *)max_min;
    }
}

/* Scaling matrix with values around overflow underflow for gerq2 */
void scale_matrix_overflow_underflow_gerq2(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char)
{

    void *max_min = NULL, *scal = NULL;
    integer max_m_n = fla_max(m, n);
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(max_m_n < 100)
        {
            calculate_scale_value(datatype, scal, max_min, 8.0, imatrix_char);
        }
        else
        {
            calculate_scale_value(datatype, scal, max_min, 10.0, imatrix_char);
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        calculate_scale_value(datatype, scal, max_min, 1.0, imatrix_char);
    }
    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);

    // Free vectors
    free_vector(max_min);
    free_vector(scal);
}

/* Scale matrix with values around overflow underflow for potrf */
void scale_matrix_overflow_underflow_potrf(integer datatype, integer m, void *A, integer lda,
                                           char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);

    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, m, lda);
        calculate_scale_value(datatype, scal, max_min, 1.0, imatrix_char);
    }
    else if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, m, lda);
        if(m < 100)
        {
            calculate_scale_value(datatype, scal, max_min, 1.1, imatrix_char);
        }
        else if(m < 500)
        {
            calculate_scale_value(datatype, scal, max_min, 1.2, imatrix_char);
        }
        else
        {
            calculate_scale_value(datatype, scal, max_min, 1.3, imatrix_char);
        }
    }
    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, m, lda, i_one);

    /* Free vectors */
    free_vector(max_min);
    free_vector(scal);
}

/* Scaling matrix with values around overflow underflow for gelqf */
void scale_matrix_underflow_overflow_gelqf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(m <= 25 && n <= 25)
        {
            calculate_scale_value(datatype, scal, max_min, 5.0, imatrix_char);
        }
        else
        {
            calculate_scale_value(datatype, scal, max_min, 9.0, imatrix_char);
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        calculate_scale_value(datatype, scal, max_min, 1.0, imatrix_char);
    }

    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);

    /* free vectors */
    free_vector(max_min);
    free_vector(scal);
}

/* Scaling matrix with values around overflow underflow for geqrf */
void scale_matrix_underflow_overflow_geqrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(m <= 25 && n <= 25)
        {
            calculate_scale_value(datatype, scal, max_min, 5.0, imatrix_char);
        }
        else
        {
            calculate_scale_value(datatype, scal, max_min, 7.0, imatrix_char);
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        calculate_scale_value(datatype, scal, max_min, 1.0, imatrix_char);
    }

    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);

    /* free vectors */
    free_vector(max_min);
    free_vector(scal);
}

/* Scaling matrix with values around overflow underflow for larfg */
void scale_matrix_underflow_overflow_larfg(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(m <= 50)
        {
            calculate_scale_value(datatype, scal, max_min, 5.0, imatrix_char);
        }
        else if (m <= 200)
        {
            calculate_scale_value(datatype, scal, max_min, 9.0, imatrix_char);
        }
        else
        {
            calculate_scale_value(datatype, scal, max_min, 25.0, imatrix_char);
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        calculate_scale_value(datatype, scal, max_min, 1.0, imatrix_char);
    }

    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);

    /* free vectors */
    free_vector(max_min);
    free_vector(scal);
}