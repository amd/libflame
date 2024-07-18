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
            flt_ratio = flt_ratio / 15.00;
            dbl_ratio = dbl_ratio / 13.00;
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
    double tuning_val = 1.0;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(max_m_n < 100)
        {
            tuning_val = 8.0;
        }
        else
        {
            tuning_val = 10.0;
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        tuning_val = 1.0;
    }
    calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);

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
    double tuning_val = 1.0;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);

    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, m, lda);
        tuning_val = 1.0;
    }
    else if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, m, lda);
        if(m < 100)
        {
            tuning_val = 1.1;
        }
        else if(m < 500)
        {
            tuning_val = 1.2;
        }
        else
        {
            tuning_val = 1.3;
        }
    }
    calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);

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
    double tuning_val = 1.0;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(m <= 25 && n <= 25)
        {
            tuning_val = 5.0;
        }
        else
        {
            tuning_val = 9.0;
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        tuning_val = 1.0;
    }
    calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);

    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);

    /* free vectors */
    free_vector(max_min);
    free_vector(scal);
}

/* Scaling matrix with values around overflow underflow for getrf */
void scale_matrix_underflow_overflow_getrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    double tuning_val = 1.0;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if((m == 1) && (n == 1))
        {
            tuning_val = 14.0;
        }
        else if((m == 1 && n < 10) || (n == 1 && m < 10))
        {
            tuning_val = 30.0;
        }
        else if((m == 1 && n < 60) || (n == 1 && m < 60))
        {
            tuning_val = 60.0;
        }
        else if((m == 1 && n < 500) || (n == 1 && m < 500))
        {
            tuning_val = 150.0;
        }
        else if(m == n && m < 30)
        {
            tuning_val = 5.0;
        }
        else if(m == n && m < 50)
        {
            tuning_val = 9.0;
        }
        else if(m == n && m < 100)
        {
            tuning_val = 13.0;
        }
        else if(m == n && m < 450)
        {
            tuning_val = 50.0;
        }
        else if(m == n && m < 800)
        {
            tuning_val = 80.0;
        }
        else if(m == n && m < 900)
        {
            tuning_val = 120.0;
        }
        else if(m == n && m < 1100)
        {
            tuning_val = 150.0;
        }
        else if(m == n)
        {
            tuning_val = 400.0;
        }
        else if(m <= 10 && n <= 10)
        {
            tuning_val = 3.5;
        }
        else if(m <= 300 && n <= 300)
        {
            tuning_val = 40.0;
        }
        else if(m <= 1000 && n <= 1000)
        {
            tuning_val = 55.0;
        }
        else
        {
            tuning_val = 80.0;
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, n, lda);
        tuning_val = 1.0;
    }
    calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);

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
            calculate_scale_value(datatype, scal, max_min, 12.0, imatrix_char);
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
        else if(m <= 200)
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

/* Scale matrix with values around overflow underflow for sytrf */
void scale_matrix_overflow_underflow_sytrf(integer datatype, integer m, void *A, integer lda,
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
            calculate_scale_value(datatype, scal, max_min, 2, imatrix_char);
        }
        else if(m < 500)
        {
            calculate_scale_value(datatype, scal, max_min, 2.1, imatrix_char);
        }
        else
        {
            calculate_scale_value(datatype, scal, max_min, 2.2, imatrix_char);
        }
    }
    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, m, lda, i_one);

    /* Free vectors */
    free_vector(max_min);
    free_vector(scal);
}

/* Scaling matrix with values around overflow underflow for org2r/ung2r */
void scale_matrix_underflow_overflow_org2r(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);
        if(m <= 50 && n <= 50)
        {
            calculate_scale_value(datatype, scal, max_min, 5.0, imatrix_char);
        }
        else if(m <= 200 && n <= 200)
        {
            calculate_scale_value(datatype, scal, max_min, 12.0, imatrix_char);
        }
        else
        {
            calculate_scale_value(datatype, scal, max_min, 35.0, imatrix_char);
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

/* Calculating the scaling value with respect to max and min for gtsv */
void calculate_gtsvA_scale_value(integer datatype, integer m, integer n, void *A, integer lda,
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
        if(fla_max(m, n) <= 25)
        {
            flt_ratio = flt_ratio / 100.0;
            dbl_ratio = dbl_ratio / 110.0;
        }
        else if(fla_max(m, n) <= 50)
        {
            flt_ratio = flt_ratio / 150.00;
            dbl_ratio = dbl_ratio / 160.00;
        }
        else if(fla_max(m, n) <= 75)
        {
            flt_ratio = flt_ratio / 200.00;
            dbl_ratio = dbl_ratio / 210.00;
        }
        else if(fla_max(m, n) <= 100)
        {
            flt_ratio = flt_ratio / 250.00;
            dbl_ratio = dbl_ratio / 260.00;
        }
        else if(fla_max(m, n) <= 150)
        {
            flt_ratio = flt_ratio / 300.00;
            dbl_ratio = dbl_ratio / 310.00;
        }
        else if(fla_max(m, n) <= 200)
        {
            flt_ratio = flt_ratio / 340.00;
            dbl_ratio = dbl_ratio / 350.00;
        }
        else if(fla_max(m, n) <= 300)
        {
            flt_ratio = flt_ratio / 430.00;
            dbl_ratio = dbl_ratio / 440.00;
        }
        else if(fla_max(m, n) <= 400)
        {
            flt_ratio = flt_ratio / 510.00;
            dbl_ratio = dbl_ratio / 520.00;
        }
        else if(fla_max(m, n) <= 500)
        {
            flt_ratio = flt_ratio / 560.00;
            dbl_ratio = dbl_ratio / 570.00;
        }
        else if(fla_max(m, n) <= 750)
        {
            flt_ratio = flt_ratio / 650.00;
            dbl_ratio = dbl_ratio / 660.00;
        }
        else if(fla_max(m, n) <= 1000)
        {
            flt_ratio = flt_ratio / 740.00;
            dbl_ratio = dbl_ratio / 750.00;
        }
        else
        {
            flt_ratio = flt_ratio / 800.00;
            dbl_ratio = dbl_ratio / 810.00;
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

/* Scaling matrix with values around overflow/underflow for gtsv */
void init_matrix_overflow_underflow_gtsv(integer datatype, integer m, integer n, void *A,
                                         integer lda, char imatrix, void *scal)
{
    calculate_gtsvA_scale_value(datatype, m, n, A, lda, imatrix, scal);
    /* Scaling the matrix A by x scalar */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);
}

/* Scaling matrix with values around overflow underflow for gels */
void scale_matrix_underflow_overflow_gels(integer datatype, char *trans, integer m, integer n,
                                          void *A, integer lda, char imatrix_char, integer sysmat)
{
    void *max_min = NULL, *scal = NULL;
    double tuning_val = 1.0;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);

    if(sysmat)
    {
        /* Tuning for A in AX=B */
        if(imatrix_char == 'O')
        {
            get_max_from_matrix(datatype, A, max_min, m, n, lda);
            if((*trans == 'T' || *trans == 'C') && m >= n)
            {
                if(datatype > DOUBLE)
                {
                    /* Tuning value for complex types */
                    tuning_val = 4.1;
                }
                else
                {
                    /* Tuning for real types */
                    tuning_val = 3.6;
                }
            }
            else
            {
                if(datatype > DOUBLE)
                {
                    /* Tuning value for complex types */
                    tuning_val = 3.4;
                }
                else
                {
                    tuning_val = 2.4;
                }
            }
        }
        if(imatrix_char == 'U')
        {
            get_min_from_matrix(datatype, A, max_min, m, n, lda);
            tuning_val = 1.0;
        }
    }
    else
    {
        /* Tuning for B in AX=B */
        if(imatrix_char == 'O')
        {
            get_max_from_matrix(datatype, A, max_min, m, n, lda);
            tuning_val = 3.4;
        }
        if(imatrix_char == 'U')
        {
            get_min_from_matrix(datatype, A, max_min, m, n, lda);
            tuning_val = 1.0;
        }
    }
    calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);

    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, n, lda, i_one);

    /* free vectors */
    free_vector(max_min);
    free_vector(scal);
}

void scale_matrix_overflow_underflow_gelsd(integer datatype, integer m, integer n, integer nrhs,
                                           void *A, integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    double tuning_val;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);

        if((m < 10) && (n < 10))
        {
            tuning_val = 25.0;
        }
        else if((m < 20) && (n < 20))
        {
            tuning_val = 40.0;
        }
        else if((m < 40) && (n < 40))
        {
            tuning_val = 48.0;
        }
        else if((m < 50) && (n < 50))
        {
            tuning_val = 55.0;
        }
        else if((m < 80) && (n < 80))
        {
            tuning_val = 70.0;
        }
        else if((m < 100) && (n < 100))
        {
            tuning_val = 78.0;
        }
        else if((m < 150) && (n < 150))
        {
            tuning_val = 98.0;
        }
        else if((m < 200) && (n < 200))
        {
            tuning_val = 120.0;
        }
        else
        {
            tuning_val = 180.0;
        }
        calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);
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

/* Scaling matrix with values around overflow underflow for gbtrf */
void scale_matrix_underflow_overflow_gbtrf(integer datatype, integer m, integer n, void *A,
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
            calculate_scale_value(datatype, scal, max_min, 4.0, imatrix_char);
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

    /* free vectors */
    free_vector(max_min);
    free_vector(scal);
}

/* Scaling matrix with values around overflow underflow for gbtrs */
void scale_matrix_underflow_overflow_gbtrs(integer datatype, integer m, integer nrhs, void *A,
                                           integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    double tuning_val = 1.0;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, nrhs, lda);
        if(m < 100)
        {
            tuning_val = 70.0;
        }
        else if(m < 600)
        {
            tuning_val = 400.0;
        }
        else if (m < 1000)
        {
            tuning_val = 700.0;
        }
        else
        {
            tuning_val = 900.0;
        }
    }
    if(imatrix_char == 'U')
    {
        get_min_from_matrix(datatype, A, max_min, m, nrhs, lda);
        tuning_val = 1.0;
    }
    calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);

    /* Scaling the matrix A with scal */
    scal_matrix(datatype, scal, A, m, nrhs, lda, i_one);

    /* free vectors */
    free_vector(max_min);
    free_vector(scal);
}

/* Scaling matrix with values around overflow underflow for gelss */
void scale_matrix_overflow_underflow_gelss(integer datatype, integer m, integer n, integer nrhs, void *A,
                                           integer lda, char imatrix_char)
{
    void *max_min = NULL, *scal = NULL;
    double tuning_val;
    create_vector(get_realtype(datatype), &max_min, 1);
    create_vector(get_realtype(datatype), &scal, 1);
    if(imatrix_char == 'O')
    {
        get_max_from_matrix(datatype, A, max_min, m, n, lda);

        if((m < 10) && (n < 10))
        {
            tuning_val = 23.0;
        }
        else if((m < 20) && (n < 20))
        {
            tuning_val = 35.0;
        }
        else if((m < 40) && (n < 40))
        {
            tuning_val = 48.0;
        }
        else if((m < 50) && (n < 50))
        {
            tuning_val = 56.0;
        }
        else if((m < 80) && (n < 80))
        {
            tuning_val = 70.0;
        }
        else if((m < 100) && (n < 100))
        {
            tuning_val = 80.0;
        }
        else if((m < 150) && (n < 150))
        {
            tuning_val = 95.0;
        }
        else if((m < 200) && (n < 200))
        {
            tuning_val = 120.0;
        }
        else
        {
            tuning_val = 180.0;
        }
        calculate_scale_value(datatype, scal, max_min, tuning_val, imatrix_char);
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

/* Scaling matrix with values around overflow, underflow for STEDC */
void scale_matrix_underflow_overflow_stedc(integer datatype, integer n, void *A, integer lda,
                                           char *imatrix_char, char *scal)
{
    float flt_ratio;
    double dbl_ratio;

    if((imatrix_char == NULL) || (scal == NULL))
    {
        return;
    }
    if(n == 1)
    {
        assign_value(datatype, A, d_one, d_one);
    }
    if(*imatrix_char == 'O')
    {
        /*Intialize the ratios with maximum of datatype value*/
        flt_ratio = FLT_MAX;
        dbl_ratio = DBL_MAX;
        void *max;
        create_vector(get_realtype(datatype), &max, 1);
        get_max_from_matrix(datatype, A, max, n, n, lda);
        /* The ratio is modified w.r.t dimension, to avoid inf in output */
        if(n <= 50)
        {
            flt_ratio = FLT_MAX / 5.0;
            dbl_ratio = DBL_MAX / 5.0;
        }
        else if(n <= 100)
        {
            flt_ratio = flt_ratio / 6.0;
            dbl_ratio = dbl_ratio / 6.0;
        }
        else if(n <= 500)
        {
            flt_ratio = flt_ratio / 12.0;
            dbl_ratio = dbl_ratio / 12.0;
        }
        else if(n <= 1000)
        {
            flt_ratio = flt_ratio / 15.0;
            dbl_ratio = dbl_ratio / 15.0;
        }
        else
        {
            flt_ratio = flt_ratio / 1000.0;
            dbl_ratio = dbl_ratio / 1000.0;
        }
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, max);
        free_vector(max);
    }
    else if(*imatrix_char == 'U')
    {
        void *min;
        create_vector(get_realtype(datatype), &min, 1);
        /* Get minimum value from matrix */
        get_min_from_matrix(datatype, A, min, n, n, lda);
        /* Intialize with minimum of datatype value */
        flt_ratio = FLT_MIN;
        dbl_ratio = DBL_MIN;
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, min);
        free_vector(min);
    }
    /* Scaling the matrix A by X scalar */
    scal_matrix(datatype, scal, A, n, n, lda, i_one);
}

/* Scaling matrix with values around overflow, underflow for STEVD */
void scale_matrix_underflow_overflow_stevd(integer datatype, integer n, void *A, integer lda,
                                           char *imatrix_char, char *scal)
{
    float flt_ratio;
    double dbl_ratio;

    if((imatrix_char == NULL) || (scal == NULL))
    {
        return;
    }
    if(n == 1)
    {
        assign_value(datatype, A, d_one, d_one);
    }
    if(*imatrix_char == 'O')
    {
        /*Intialize the ratios with maximum of datatype value*/
        flt_ratio = FLT_MAX;
        dbl_ratio = DBL_MAX;
        void *max;
        create_vector(get_realtype(datatype), &max, 1);
        get_max_from_matrix(datatype, A, max, n, n, lda);
        /* The ratio is modified w.r.t dimension, to avoid inf in output */
        if(n <= 50)
        {
            flt_ratio = FLT_MAX / 4.0;
            dbl_ratio = DBL_MAX / 4.0;
        }
        else if(n <= 100)
        {
            flt_ratio = flt_ratio / 4.0;
            dbl_ratio = dbl_ratio / 4.5;
        }
        else if(n <= 500)
        {
            flt_ratio = flt_ratio / 9.0;
            dbl_ratio = dbl_ratio / 9.0;
        }
        else if(n <= 1000)
        {
            flt_ratio = flt_ratio / 12.0;
            dbl_ratio = dbl_ratio / 12.0;
        }
        else
        {
            flt_ratio = flt_ratio / 1000.0;
            dbl_ratio = dbl_ratio / 1000.0;
        }
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, max);
        free_vector(max);
    }
    else if(*imatrix_char == 'U')
    {
        void *min;
        create_vector(get_realtype(datatype), &min, 1);
        /* Get minimum value from matrix */
        get_min_from_matrix(datatype, A, min, n, n, lda);
        /* Intialize with minimum of datatype value */
        flt_ratio = FLT_MIN;
        dbl_ratio = DBL_MIN;
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, min);
        free_vector(min);
    }
    /* Scaling the matrix A by X scalar */
    scal_matrix(datatype, scal, A, n, n, lda, i_one);
}

/* Scaling matrix with values around overflow, underflow for SYEV/HEEV */
void scale_matrix_underflow_overflow_syev(integer datatype, integer n, void *A, integer lda,
                                          char *imatrix_char, char *scal)
{
    float flt_ratio;
    double dbl_ratio;

    if((imatrix_char == NULL) || (scal == NULL))
    {
        return;
    }
    if(n == 1)
    {
        assign_value(datatype, A, d_one, d_one);
    }
    if(*imatrix_char == 'O')
    {
        /*Intialize the ratios with maximum of datatype value*/
        flt_ratio = FLT_MAX;
        dbl_ratio = DBL_MAX;
        void *max;
        create_vector(get_realtype(datatype), &max, 1);
        get_max_from_matrix(datatype, A, max, n, n, lda);
        /* The ratio is modified w.r.t dimension, to avoid inf in output */
        if(n <= 50)
        {
            flt_ratio = FLT_MAX / 4.5;
            dbl_ratio = DBL_MAX / 5.0;
        }
        else if(n <= 100)
        {
            flt_ratio = flt_ratio / 6.0;
            dbl_ratio = dbl_ratio / 6.0;
        }
        else if(n <= 500)
        {
            flt_ratio = flt_ratio / 12.0;
            dbl_ratio = dbl_ratio / 12.0;
        }
        else if(n <= 1000)
        {
            flt_ratio = flt_ratio / 16.0;
            dbl_ratio = dbl_ratio / 16.0;
        }
        else
        {
            flt_ratio = flt_ratio / 1000.0;
            dbl_ratio = dbl_ratio / 1000.0;
        }
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, max);
        free_vector(max);
    }
    else if(*imatrix_char == 'U')
    {
        void *min;
        create_vector(get_realtype(datatype), &min, 1);
        /* Get minimum value from matrix */
        get_min_from_matrix(datatype, A, min, n, n, lda);
        /* Intialize with minimum of datatype value */
        flt_ratio = FLT_MIN;
        dbl_ratio = DBL_MIN;
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, min);
        free_vector(min);
    }
    /* Scaling the matrix A by X scalar */
    scal_matrix(datatype, scal, A, n, n, lda, i_one);
}

/* Scaling matrix with values around overflow, underflow for STEQR */
void scale_matrix_underflow_overflow_steqr(integer datatype, integer n, void *A, integer lda,
                                           char *imatrix_char, char *scal)
{
    float flt_ratio;
    double dbl_ratio;

    if((imatrix_char == NULL) || (scal == NULL))
    {
        return;
    }
    if(n == 1)
    {
        assign_value(datatype, A, d_one, d_one);
    }
    if(*imatrix_char == 'O')
    {
        /*Intialize the ratios with maximum of datatype value*/
        flt_ratio = FLT_MAX;
        dbl_ratio = DBL_MAX;
        void *max;
        create_vector(get_realtype(datatype), &max, 1);
        get_max_from_matrix(datatype, A, max, n, n, lda);
        /* The ratio is modified w.r.t dimension, to avoid inf in output */
        if(n <= 50)
        {
            flt_ratio = FLT_MAX / 6.0;
            dbl_ratio = DBL_MAX / 6.0;
        }
        else if(n <= 100)
        {
            flt_ratio = flt_ratio / 7.0;
            dbl_ratio = dbl_ratio / 7.0;
        }
        else if(n <= 500)
        {
            flt_ratio = flt_ratio / 13.0;
            dbl_ratio = dbl_ratio / 13.5;
        }
        else if(n <= 1000)
        {
            flt_ratio = flt_ratio / 16.0;
            dbl_ratio = dbl_ratio / 17.0;
        }
        else
        {
            flt_ratio = flt_ratio / 1000.0;
            dbl_ratio = dbl_ratio / 1000.0;
        }
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, max);
        free_vector(max);
    }
    else if(*imatrix_char == 'U')
    {
        void *min;
        create_vector(get_realtype(datatype), &min, 1);
        /* Get minimum value from matrix */
        get_min_from_matrix(datatype, A, min, n, n, lda);
        /* Intialize with minimum of datatype value */
        flt_ratio = FLT_MIN;
        dbl_ratio = DBL_MIN;
        compute_ratio(get_realtype(datatype), scal, flt_ratio, dbl_ratio, min);
        free_vector(min);
    }
    /* Scaling the matrix A by X scalar */
    scal_matrix(datatype, scal, A, n, n, lda, i_one);
}
