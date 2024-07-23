/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_larf.c
 *  @brief Defines validate function of LARF() to use in test suite.
 *  */
#include "test_common.h"

void validate_larf(integer datatype, char side, integer m, integer n, void *v, integer incv,
                   void *c__, integer ldc__, void *c__out, integer ldc__out, void *tau,
                   double *residual)
{
    /*
    This function does the following:
    1.1 If type = float or double, compute H = I - tau * v_tmp * v_tmp-transpose
    1.2 If type = complex or double complex, compute H = I - tau * v_tmp * v_tmp-hermitian-transpose

    2.1 If type = float or double,  compute out_validate = H * c__out_tmp(when side = L) or
    c_out_tmp * H(when side = R)
    2.2 If type = complex or double complex,  compute out_validate =
    H-hermitian-transpose * c__out_tmp(when side = L) or c_out_tmp * H-hermitian-transpose(when side
    = R)

    3. Compute residual using c__tmp and out_validate
    */

    void *v_tmp;
    void *H;
    void *c__out_tmp;
    void *out_validate;
    void *work;
    void *c__tmp;

    int m_H;
    if(side == 'L')
    {
        m_H = m;
    }
    else
    {
        m_H = n;
    }

    create_vector(datatype, &v_tmp, m_H);
    create_matrix(datatype, matrix_layout, m_H, m_H, &H, m_H);
    set_identity_matrix(datatype, m_H, m_H, H, m_H);
    copy_vector(datatype, m_H, v, incv, v_tmp, 1);

    switch(datatype)
    {
        case FLOAT:
        {
            /* 1. Compute H = I - tau * v_tmp * v_tmp-transpose */
            *((float *)tau) = -(*((float *)tau));
            sger_(&m_H, &m_H, tau, v_tmp, &i_one, v_tmp, &i_one, H, &m_H);

            /* 2. Compute out_validate = H * c__out_tmp(side = L) or c_out_tmp * H(side = R)*/
            create_matrix(datatype, matrix_layout, m, n, &c__out_tmp, m);
            create_matrix(datatype, matrix_layout, m, n, &out_validate, m);

            copy_matrix(datatype, "full", m, n, c__out, ldc__out, c__out_tmp, m);

            if(side == 'L')
            {
                sgemm_("N", "N", &m, &n, &m, &s_one, H, &m, c__out_tmp, &m, &s_zero, out_validate,
                       &m);
            }
            else
            {
                sgemm_("N", "N", &m, &n, &n, &s_one, c__out_tmp, &m, H, &n, &s_zero, out_validate,
                       &m);
            }

            /* 3. Compute residual using c__tmp and out_validate */
            create_matrix(datatype, matrix_layout, m, n, &c__tmp, m);

            copy_matrix(datatype, "full", m, n, c__, ldc__, c__tmp, m);
            create_vector(datatype, &work, m);

            float norm, norm1, eps, resid1;

            eps = fla_lapack_slamch("P");
            norm1 = fla_lapack_slange("1", &m, &n, c__tmp, &m, work);
            matrix_difference(datatype, m, n, c__tmp, m, out_validate, m);
            norm = fla_lapack_slange("1", &m, &n, c__tmp, &m, work);

            resid1 = norm / (float)m / norm1 / eps;
            *residual = (double)resid1;
            break;
        }
        case DOUBLE:
        {
            /* 1. Compute H = I - tau * v_tmp * v_tmp-transpose */
            *((double *)tau) = -(*((double *)tau));
            dger_(&m_H, &m_H, tau, v_tmp, &i_one, v_tmp, &i_one, H, &m_H);

            /* 2. Compute out_validate = H * c__out_tmp(side = L) or c_out_tmp * H(side = R) */
            create_matrix(datatype, matrix_layout, m, n, &c__out_tmp, m);
            create_matrix(datatype, matrix_layout, m, n, &out_validate, m);

            copy_matrix(datatype, "full", m, n, c__out, ldc__out, c__out_tmp, m);

            if(side == 'L')
            {
                dgemm_("N", "N", &m, &n, &m, &d_one, H, &m, c__out_tmp, &m, &d_zero, out_validate,
                       &m);
            }
            else
            {
                dgemm_("N", "N", &m, &n, &n, &d_one, c__out_tmp, &m, H, &n, &d_zero, out_validate,
                       &m);
            }

            /* 3. Compute residual using c__tmp and out_validate */
            create_matrix(datatype, matrix_layout, m, n, &c__tmp, m);

            copy_matrix(datatype, "full", m, n, c__, ldc__, c__tmp, m);
            create_vector(datatype, &work, m);

            double norm, norm1, eps, resid1;

            eps = fla_lapack_dlamch("P");
            norm1 = fla_lapack_dlange("1", &m, &n, c__tmp, &m, work);
            matrix_difference(datatype, m, n, c__tmp, m, out_validate, m);
            norm = fla_lapack_dlange("1", &m, &n, c__tmp, &m, work);

            resid1 = norm / (double)m / norm1 / eps;
            *residual = resid1;
            break;
        }
        case COMPLEX:
        {
            /* 1. Compute H = I - tau * v_tmp * v_tmp-hermitian-transpose */
            ((scomplex *)tau)->real = -(((scomplex *)tau)->real);
            ((scomplex *)tau)->imag = -(((scomplex *)tau)->imag);
            cgerc_(&m_H, &m_H, tau, v_tmp, &i_one, v_tmp, &i_one, H, &m_H);

            /* 2. Compute out_validate = H-hermitian-transpose * c__out_tmp(when side = L) or
             * c_out_tmp * H-hermitian-transpose(when side = R) */
            create_matrix(datatype, matrix_layout, m, n, &c__out_tmp, m);
            create_matrix(datatype, matrix_layout, m, n, &out_validate, m);

            copy_matrix(datatype, "full", m, n, c__out, ldc__out, c__out_tmp, m);

            if(side == 'L')
            {
                cgemm_("C", "N", &m, &n, &m, &c_one, H, &m, c__out_tmp, &m, &c_zero, out_validate,
                       &m);
            }
            else
            {
                cgemm_("N", "C", &m, &n, &n, &c_one, c__out_tmp, &m, H, &n, &c_zero, out_validate,
                       &m);
            }

            /* 3. Compute residual using c__tmp and out_validate */
            create_matrix(datatype, matrix_layout, m, n, &c__tmp, m);

            copy_matrix(datatype, "full", m, n, c__, ldc__, c__tmp, m);
            create_vector(datatype, &work, m);

            float norm, norm1, eps, resid1;

            eps = fla_lapack_slamch("P");
            norm1 = fla_lapack_clange("1", &m, &n, c__tmp, &m, work);
            matrix_difference(datatype, m, n, c__tmp, m, out_validate, m);
            norm = fla_lapack_clange("1", &m, &n, c__tmp, &m, work);

            resid1 = norm / (float)m / norm1 / eps;
            *residual = (double)resid1;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* 1. Compute H = I - tau * v_tmp * v_tmp-hermitian-transpose */
            ((dcomplex *)tau)->real = -(((dcomplex *)tau)->real);
            ((dcomplex *)tau)->imag = -(((dcomplex *)tau)->imag);
            zgerc_(&m_H, &m_H, tau, v_tmp, &i_one, v_tmp, &i_one, H, &m_H);

            /* 2. Compute out_validate = H-hermitian-transpose * c__out_tmp(when side = L) or
             * c_out_tmp * H-hermitian-transpose(when side = R) */
            create_matrix(datatype, matrix_layout, m, n, &c__out_tmp, m);
            create_matrix(datatype, matrix_layout, m, n, &out_validate, m);

            copy_matrix(datatype, "full", m, n, c__out, ldc__out, c__out_tmp, m);

            if(side == 'L')
            {
                zgemm_("C", "N", &m, &n, &m, &z_one, H, &m, c__out_tmp, &m, &z_zero, out_validate,
                       &m);
            }
            else
            {
                zgemm_("N", "C", &m, &n, &n, &z_one, c__out_tmp, &m, H, &n, &z_zero, out_validate,
                       &m);
            }

            /* 3. Compute residual using c__tmp and out_validate */
            create_matrix(datatype, matrix_layout, m, n, &c__tmp, m);

            copy_matrix(datatype, "full", m, n, c__, ldc__, c__tmp, m);
            create_vector(datatype, &work, m);

            double norm, norm1, eps, resid1;

            eps = fla_lapack_dlamch("P");
            norm1 = fla_lapack_zlange("1", &m, &n, c__tmp, &m, work);
            matrix_difference(datatype, m, n, c__tmp, m, out_validate, m);
            norm = fla_lapack_zlange("1", &m, &n, c__tmp, &m, work);

            resid1 = norm / (double)m / norm1 / eps;
            *residual = resid1;
            break;
        }
    }

    /* Free the matrices and vectors */

    free_vector(v_tmp);
    free_matrix(H);
    free_matrix(c__out_tmp);
    free_matrix(out_validate);
    free_matrix(c__tmp);
    free_vector(work);
}
