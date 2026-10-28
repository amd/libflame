/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gebal.c
 *  @brief Defines validate function of GEBAL() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

#include <math.h>

extern double perf;
extern double time_min;

/* JOB N and JOB P must not scale the active block. */
static integer gebal_active_scale_is_one(integer datatype, integer ilo, integer ihi, void *scale)
{
    integer scale_datatype = get_realtype(datatype);
    integer i;

    for(i = ilo - 1; i < ihi; i++)
    {
        if(get_realtype_value(scale_datatype, get_ptr_at_offset(scale_datatype, scale, i)) != 1.0)
            return 0;
    }
    return 1;
}

/* JOB S and JOB B must return finite positive scale factors for the active block. */
static integer gebal_active_scale_is_valid(integer datatype, integer ilo, integer ihi, void *scale)
{
    integer scale_datatype = get_realtype(datatype);
    integer i;

    for(i = ilo - 1; i < ihi; i++)
    {
        double value = get_realtype_value(scale_datatype, get_ptr_at_offset(scale_datatype, scale, i));

        if(!isfinite(value) || value <= 0.0)
            return 0;
    }
    return 1;
}

/* SCALE holds a 1-based permutation index outside [ILO, IHI]. Reject
   non-finite, non-integral, and out-of-range values before converting. */
static integer gebal_partner(double value, integer n, integer *ok)
{
    integer p;

    if(!isfinite(value) || value < 1.0 || value > (double)n || floor(value) != value)
    {
        *ok = 0;
        return 0;
    }
    p = (integer)value;
    if(p < 1 || p > n)
    {
        *ok = 0;
        return 0;
    }
    return p;
}

/* Return x * r * 2^e without forming the full ratio. The exponent is applied in two
   halves around the mantissa multiply, so an overflowing di/dj (e.g. long double equal
   to double on Microsoft-ABI builds) is never materialized. */
static double gebal_scale_value(double x, double r, int e)
{
    int e1 = e / 2;

    return scalbn(scalbn(x, e1) * r, e - e1);
}

/* Multiply one matrix entry by (rm * 2^re), a real balancing ratio kept as a mantissa
   and a power-of-two exponent. */
static void gebal_mul(integer datatype, void *A, integer lda, integer i, integer j, double rm,
                      int re)
{
    switch(datatype)
    {
        case FLOAT:
            ((float *)A)[j * lda + i]
                = (float)gebal_scale_value((double)((float *)A)[j * lda + i], rm, re);
            break;
        case DOUBLE:
            ((double *)A)[j * lda + i] = gebal_scale_value(((double *)A)[j * lda + i], rm, re);
            break;
        case COMPLEX:
            ((scomplex *)A)[j * lda + i].real
                = (float)gebal_scale_value((double)((scomplex *)A)[j * lda + i].real, rm, re);
            ((scomplex *)A)[j * lda + i].imag
                = (float)gebal_scale_value((double)((scomplex *)A)[j * lda + i].imag, rm, re);
            break;
        case DOUBLE_COMPLEX:
            ((dcomplex *)A)[j * lda + i].real
                = gebal_scale_value(((dcomplex *)A)[j * lda + i].real, rm, re);
            ((dcomplex *)A)[j * lda + i].imag
                = gebal_scale_value(((dcomplex *)A)[j * lda + i].imag, rm, re);
            break;
    }
}

/* Undo diagonal scaling. Balanced(i,j) = Permuted(i,j) * D(j) / D(i),
   with D = 1 outside [ilo, ihi]. Each D is split with frexp into mantissa * 2^exp, so
   the ratio is mantissa ratio (in [0.25, 4]) and an exponent difference. */
static integer gebal_undo_scale(integer datatype, integer n, integer ilo, integer ihi, void *R,
                                integer lda, void *scale)
{
    integer scale_datatype = get_realtype(datatype);
    integer i, j;

    for(j = 0; j < n; j++)
    {
        integer jb = j + 1;
        double dj = 1.0, mj;
        int ej;

        if(jb >= ilo && jb <= ihi)
            dj = get_realtype_value(scale_datatype, get_ptr_at_offset(scale_datatype, scale, j));
        if(!isfinite(dj) || dj == 0.0)
            return 0;
        mj = frexp(dj, &ej);

        for(i = 0; i < n; i++)
        {
            integer ib = i + 1;
            double di = 1.0, mi;
            int ei;

            if(ib >= ilo && ib <= ihi)
                di = get_realtype_value(scale_datatype,
                                        get_ptr_at_offset(scale_datatype, scale, i));
            if(!isfinite(di) || di == 0.0)
                return 0;
            mi = frexp(di, &ei);
            gebal_mul(datatype, R, lda, i, j, mi / mj, ei - ej);
        }
    }
    return 1;
}

/* Permutations were applied from N down to IHI+1, then from 1 up to ILO-1.
   Each step swapped a column, then a row. Undo in reverse, row then column. */
static integer gebal_undo_perm(integer datatype, integer n, integer ilo, integer ihi, void *R,
                               integer lda, void *scale)
{
    integer scale_datatype = get_realtype(datatype);
    integer k, j, ok = 1;
    integer cnt, inc_one = 1, inc_lda = lda;

    for(k = ilo - 1; k >= 1 && ok; k--)
    {
        integer p = gebal_partner(
            get_realtype_value(scale_datatype, get_ptr_at_offset(scale_datatype, scale, k - 1)), n,
            &ok);
        if(!ok || p == k)
            continue;
        cnt = n - (k - 1);
        swap_row_col(datatype, &cnt, R, lda, &inc_lda, &inc_lda, p - 1, k - 1, k - 1, k - 1);
        swap_row_col(datatype, &ihi, R, lda, &inc_one, &inc_one, 0, p - 1, 0, k - 1);
    }
    for(j = ihi + 1; j <= n && ok; j++)
    {
        integer p = gebal_partner(
            get_realtype_value(scale_datatype, get_ptr_at_offset(scale_datatype, scale, j - 1)), n,
            &ok);
        if(!ok || p == j)
            continue;
        swap_row_col(datatype, &n, R, lda, &inc_lda, &inc_lda, p - 1, 0, j - 1, 0);
        swap_row_col(datatype, &j, R, lda, &inc_one, &inc_one, 0, p - 1, 0, j - 1);
    }
    return ok;
}

/* Permutation bounds derived only from A. This does not read ILO, IHI, or
   SCALE from the routine under test, so a no-op that reports the full
   interval fails when A has a row or column that must be isolated. */
static void gebal_expected_perm_bounds(integer datatype, integer n, void *A, integer lda,
                                       integer *exp_ilo, integer *exp_ihi)
{
    void *B = NULL;
    integer i, j, k, l;
    integer noconv, canswap;
    integer cnt, inc_one = 1, inc_lda = lda;

    *exp_ilo = 1;
    *exp_ihi = n;
    if(n <= 0)
        return;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B, lda);
    copy_matrix(datatype, "full", n, n, A, lda, B, lda);

    k = 1;
    l = n;
    noconv = 1;
    while(noconv)
    {
        noconv = 0;
        for(i = l; i >= 1; --i)
        {
            canswap = 1;
            for(j = 1; j <= l; ++j)
            {
                if(i != j && get_abs(datatype, B, lda, i - 1, j - 1) != 0.0)
                {
                    canswap = 0;
                    break;
                }
            }
            if(!canswap)
                continue;
            if(i != l)
            {
                swap_row_col(datatype, &l, B, lda, &inc_one, &inc_one, 0, i - 1, 0, l - 1);
                cnt = n - (k - 1);
                swap_row_col(datatype, &cnt, B, lda, &inc_lda, &inc_lda, i - 1, k - 1, l - 1,
                             k - 1);
            }
            noconv = 1;
            if(l == 1)
            {
                *exp_ilo = 1;
                *exp_ihi = 1;
                free_matrix(B);
                return;
            }
            --l;
        }
    }

    noconv = 1;
    while(noconv)
    {
        noconv = 0;
        for(j = k; j <= l; ++j)
        {
            canswap = 1;
            for(i = k; i <= l; ++i)
            {
                if(i != j && get_abs(datatype, B, lda, i - 1, j - 1) != 0.0)
                {
                    canswap = 0;
                    break;
                }
            }
            if(!canswap)
                continue;
            if(j != k)
            {
                cnt = l;
                swap_row_col(datatype, &cnt, B, lda, &inc_one, &inc_one, 0, j - 1, 0, k - 1);
                cnt = n - (k - 1);
                swap_row_col(datatype, &cnt, B, lda, &inc_lda, &inc_lda, j - 1, k - 1, k - 1,
                             k - 1);
            }
            noconv = 1;
            ++k;
        }
    }

    *exp_ilo = k;
    *exp_ihi = l;
    free_matrix(B);
}

/* 1-norm of entries that GEBAL must leave zero outside the balanced block. */
static double gebal_structure_norm(integer datatype, integer n, integer ilo, integer ihi, void *A,
                                   integer lda)
{
    integer i, j;
    double struct_norm = 0.0;

    for(j = 0; j < n; j++)
    {
        integer jb = j + 1;
        double colsum = 0.0;

        for(i = 0; i < n; i++)
        {
            integer ib = i + 1;
            integer illegal = 0;

            if(jb < ilo && ib > jb)
                illegal = 1;
            if(ib > ihi && jb < ib)
                illegal = 1;
            if(illegal)
                colsum += get_abs(datatype, A, lda, i, j);
        }
        if(colsum > struct_norm)
            struct_norm = colsum;
    }
    return struct_norm;
}

/* Validate GEBAL by checking JOB semantics, permutation structure, inverse
   scaling/permutation reconstruction, and padding preservation. */
void validate_gebal(char *tst_api, char job, integer n, integer ilo, integer ihi, void *A,
                    void *A_test, integer lda, void *scale, integer datatype,
                    double err_thresh, void *params)
{
    void *R = NULL, *work = NULL;
    integer recons_check, perm_check = 1, bounds_valid, semantics_valid = 1;
    double residual, resid1 = 0., resid2 = 0., resid3 = 0.;
    double norm = 0., norm_A = 0.;
    float nrm_s = 0.0f;
    double nrm_d = 0.0;
    integer single_prec = (datatype == FLOAT || datatype == COMPLEX);
    void *nrm = single_prec ? (void *)&nrm_s : (void *)&nrm_d;

    if(n == 0)
    {
        /* For n == 0, the routine should set an empty active block (ILO=1, IHI=0). */
        if(!(ilo == 1 && ihi == 0))
        {
            residual = DBL_MAX;
            FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
            return;
        }
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &R, lda);
    create_realtype_vector(datatype, &work, n);

    /* 1-norm of the original matrix, used to normalize both residuals. */
    compute_matrix_norm(datatype, '1', n, n, A, lda, nrm, 'N', work);
    norm_A = single_prec ? (double)nrm_s : nrm_d;

    /* ILO/IHI must address SCALE. An empty active block is ILO = IHI + 1.
       JOB 'N' and 'S' do not permute. JOB 'P' and 'B' must match the
       permutation interval computed from the original matrix. */
    bounds_valid = (ilo >= 1 && ihi >= 1 && ihi <= n && ilo <= ihi + 1);
    if(!bounds_valid)
        perm_check = 0;
    else if((same_char(job, 'N') || same_char(job, 'S')) && (ilo != 1 || ihi != n))
        perm_check = 0;
    else if(same_char(job, 'P') || same_char(job, 'B'))
    {
        integer exp_ilo, exp_ihi;

        gebal_expected_perm_bounds(datatype, n, A, lda, &exp_ilo, &exp_ihi);
        if(ilo != exp_ilo || ihi != exp_ihi)
            perm_check = 0;
    }

    /* Enforce the operation selected by JOB, rather than accepting any
       balancing transformation that happens to reconstruct A. */
    if(!bounds_valid)
        semantics_valid = 0;
    else if(same_char(job, 'N') || same_char(job, 'P'))
        semantics_valid = gebal_active_scale_is_one(datatype, ilo, ihi, scale);
    else if(same_char(job, 'S') || same_char(job, 'B'))
    {
        semantics_valid = gebal_active_scale_is_valid(datatype, ilo, ihi, scale);
    }
    else
        semantics_valid = 0;

    if(!perm_check || !semantics_valid)
        resid2 = DBL_MAX;
    else
        resid2 = fla_compute_residual(datatype, 'P',
                                      gebal_structure_norm(datatype, n, ilo, ihi, A_test, lda),
                                      norm_A, n, params);

    copy_matrix(datatype, "full", n, n, A_test, lda, R, lda);

    if(!bounds_valid)
        recons_check = 0;
    else
    {
        recons_check = gebal_undo_scale(datatype, n, ilo, ihi, R, lda, scale);
        if(recons_check)
            recons_check = gebal_undo_perm(datatype, n, ilo, ihi, R, lda, scale);
    }

    if(!recons_check)
    {
        resid1 = DBL_MAX;
    }
    else
    {
        matrix_difference(datatype, n, n, R, lda, A, lda);
        compute_matrix_norm(datatype, '1', n, n, R, lda, nrm, 'N', work);
        norm = single_prec ? (double)nrm_s : nrm_d;
        if(norm_A == 0.0)
            resid1 = (norm == 0.0) ? 0.0 : DBL_MAX;
        else
        {
            resid1 = fla_compute_residual(datatype, 'P', norm, norm_A, n, params);
            /* If the helper returned zero due to a tiny base but the difference is
               nonzero, fail explicitly instead of silently accepting the output. */
            if(resid1 == 0.0 && norm != 0.0)
                resid1 = DBL_MAX;
        }
    }

    resid3 = check_padding(datatype, n, n, A_test, lda);

    free_vector(work);
    free_matrix(R);

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(residual, resid3);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
}
