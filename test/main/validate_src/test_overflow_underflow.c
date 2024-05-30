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

/* Finding ratio between two real datatype values */
void compute_ratio(integer datatype, void *scal, float flt_quotient, double dbl_quotient,
                   void *denominator)
{
    if(get_realtype(datatype) == FLOAT)
        *(float *)scal = flt_quotient / *(float *)denominator;
    else if(get_realtype(datatype) == DOUBLE)
        *(double *)scal = dbl_quotient / *(double *)denominator;
}

/* Initializing matrix with values around overflow underflow */
void init_matrix_overflow_underflow_svd(integer datatype, integer m, integer n, void *A,
                                        integer lda, char imatrix, void *scal)
{
    /* Find factor to scale up/down the matrix */
    calculate_svd_scale_value(datatype, m, n, A, lda, imatrix, scal);
    /* Scaling the matrix A by x scalar */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);
}

/* Scaling matrix with values around overflow underflow for gerq2 */
void scale_matrix_overflow_underflow_gerq2(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char)
{
    float feasible_float = 0.0;
    double feasible_double = 0.0;
    void *max_min = NULL, *scal = NULL;
    integer max_m_n = fla_max(m, n);
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        feasible_float = FLT_MAX;
        feasible_double = DBL_MAX;
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(max_m_n < 100)
        {
            if(get_realtype(datatype) == FLOAT)
            {
                *(float *)scal = (feasible_float / 8.0 ) / *(float *)max_min;
            }
            else if(get_realtype(datatype) == DOUBLE)
            {
                *(double *)scal = (feasible_double / 8.0) / *(double *)max_min;
            }
        }
        else
        {
            if(get_realtype(datatype) == FLOAT)
            {
                *(float *)scal = (feasible_float / 10.0) / *(float *)max_min;
            }
            else if(get_realtype(datatype) == DOUBLE)
            {
                *(double *)scal = (feasible_double / 10.0) / *(double *)max_min;
            }
        }
    }
    if(imatrix_char == 'U')
    {
        feasible_float = FLT_MIN;
        feasible_double = DBL_MIN;
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        if(get_realtype(datatype) == FLOAT)
        {
            *(float *)scal = feasible_float / *(float *)max_min;
        }
        else if(get_realtype(datatype) == DOUBLE)
        {
            *(double *)scal = feasible_double / *(double *)max_min;
        }
    }
    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);

    // Free vectors
    free_vector(max_min);
    free_vector(scal);
}

/* Scale matrix with values around overflow underflow for potrf */
void scale_matrix_overflow_underflow_potrf(integer datatype, integer m, void *A, integer lda,
                                           char imatrix)
{
    float feasible_float = 0.0;
    double feasible_double = 0.0;
    void *max_min = NULL, *scal = NULL;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);

    if(imatrix == 'U')
    {
        feasible_float = FLT_MIN;
        feasible_double = DBL_MIN;
        get_min_from_matrix(datatype, A, max_min, m, m, lda);
        {
            if(get_realtype(datatype) == FLOAT)
            {
                *(float *)scal = feasible_float / *(float *)max_min;
            }
            else if(get_realtype(datatype) == DOUBLE)
            {
                *(double *)scal = feasible_double / *(double *)max_min;
            }
        }
    }
    else if(imatrix == 'O')
    {
        feasible_float = FLT_MAX;
        feasible_double = DBL_MAX;
        get_max_from_matrix(datatype, A, max_min, m, m, lda);
        if(m < 100)
        {
            if(get_realtype(datatype) == FLOAT)
            {
                *(float *)scal = (feasible_float / *(float *)max_min) / 1.1;
            }
            else if(get_realtype(datatype) == DOUBLE)
            {
                *(double *)scal = (feasible_double / *(double *)max_min) / 1.1;
            }
        }
        else if(m < 500)
        {
            if(get_realtype(datatype) == FLOAT)
            {
                *(float *)scal = (feasible_float / *(float *)max_min) / 1.2;
            }
            else if(get_realtype(datatype) == DOUBLE)
            {
                *(double *)scal = (feasible_double / *(double *)max_min) / 1.2;
            }
        }
        else
        {
            if(get_realtype(datatype) == FLOAT)
            {
                *(float *)scal = (feasible_float / *(float *)max_min) / 1.3;
            }
            else if(get_realtype(datatype) == DOUBLE)
            {
                *(double *)scal = (feasible_double / *(double *)max_min) / 1.3;
            }
        }
    }
    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, m, lda, i_one);

    /* Free vectors */
    free_vector(max_min);
    free_vector(scal);
}