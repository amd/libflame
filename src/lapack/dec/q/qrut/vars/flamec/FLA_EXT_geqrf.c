/*
    Copyright (c) 2021-2023 Advanced Micro Devices, Inc.  All rights reserved.
    May 09, 2021
*/

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

#define ssign( x ) ( (x) < 0.0F ? -1.0F : 1.0F )
#define dsign( x ) ( (x) < 0.0  ? -1.0  : 1.0  )


FLA_Error FLA_EXT_Househ2_l_ops( fla_dim_t  m_x2,
                                 float*   chi_1,
                                 float*   x2, fla_dim_t inc_x2,
                                 float*   tau )
{
  float   one_half = 1.0F/2.0F;
  float   y[2];
  float   alpha;
  float   chi_1_minus_alpha, inv_chi_1_minus_alpha;
  float   norm_x_2;
  float   norm_x;
  float   safmin, rsafmn, lchi1;
  fla_dim_t i_one = 1;
  fla_dim_t i_two = 2;
  fla_dim_t kn;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  norm_x_2 = aocl_blas_snrm2( &m_x2,
                     x2, &inc_x2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0 )
  {
    *chi_1 = -(*chi_1);
    *tau   = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  lchi1 = *chi_1;
  y[0] = lchi1;
  y[1] = norm_x_2;

  norm_x = aocl_blas_snrm2( &i_two,
                   y, &i_one );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //          = -sign( chi_1 ) * || x ||_2
  //

  alpha = -ssign( lchi1 ) * norm_x;

  //
  // Overwrite x_2 with u_2:
  //
  //   x_2 := x_2 / ( chi_1 - alpha )
  //

  chi_1_minus_alpha = lchi1 - alpha;

  //
  // Scale X2, chi_1_minus_alpha, alpha
  // if norm factor is very less
  //
  safmin = fla_slamch("S", 1) / fla_slamch("E", 1);
  kn = 0;
  if( fabs( chi_1_minus_alpha ) < safmin )
  {
    rsafmn = 1. / safmin;

    for( kn = 1; kn < 20; kn++ )
    {
      aocl_blas_sscal( &m_x2,
              &rsafmn,
              x2, &inc_x2 );
      alpha = alpha * rsafmn;
      chi_1_minus_alpha = chi_1_minus_alpha * rsafmn;
      
      if( fabs( chi_1_minus_alpha ) > safmin )
        break;
    }
  }

  inv_chi_1_minus_alpha = 1.0F / chi_1_minus_alpha;

  aocl_blas_sscal( &m_x2,
          &inv_chi_1_minus_alpha,
          x2, &inc_x2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //        = alpha / ( alpha - chi_1 )
  //

  *tau = -1.0F * alpha * inv_chi_1_minus_alpha;

  //
  // Scale back alpha
  //
  for( ; kn > 0; kn-- )
  {
    alpha = alpha * safmin;
  }

  //
  // Overwrite chi_1 with alpha:
  //
  //   chi_1 := alpha
  //

  *chi_1 = alpha;

  return FLA_SUCCESS;
}

FLA_Error FLA_EXT_sgeqrf( fla_dim_t  m_A, fla_dim_t n_A,
                          float*   buff_A, fla_dim_t cs_A,
                          float*   buff_t,
                          float*   buff_w,
                          fla_dim_t* lwork,
                          fla_dim_t* info )
{
  fla_dim_t min_m_n = fla_min( m_A, n_A );
  fla_dim_t i, rs_A = 1;

  for ( i = 0; i < min_m_n; ++i )
  {
    float* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    float* tau1     = buff_t + i;
    float  alphat = *alpha11;

    fla_dim_t m_curr   = m_A - i;

    fla_dim_t m_ahead  = m_A - i - 1;
    fla_dim_t n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // Compute Householder transformation for current column
    FLA_EXT_Househ2_l_ops( m_ahead,
                           &alphat,
                           a21, rs_A,
                           tau1 );

    *alpha11 = 1.0F;
    if( *tau1 != 0.0F )
        *tau1 = 1.0F / *tau1;

    // Apply the computed Householder transformation on the matrix
    aocl_lapack_slarf( "Left",
            &m_curr, &n_ahead,
            alpha11, &rs_A,
            tau1,
            a12t, &cs_A,
            buff_w );

    *alpha11 = alphat;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

FLA_Error FLA_EXT_Househ2_l_opd( fla_dim_t   m_x2,
                                 double*   chi_1,
                                 double*   x2, fla_dim_t inc_x2,
                                 double*   tau )
{
  double   one_half = 1.0/2.0;
  double   y[2];
  double   alpha;
  double   chi_1_minus_alpha, inv_chi_1_minus_alpha;
  double   norm_x_2;
  double   norm_x;
  double   safmin, rsafmn, lchi1;
  fla_dim_t  i_one = 1;
  fla_dim_t  i_two = 2;
  fla_dim_t  kn;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  norm_x_2 = aocl_blas_dnrm2( &m_x2,
                     x2, &inc_x2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0 )
  {
    *chi_1 = -(*chi_1);
    *tau   = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  lchi1 = *chi_1;
  y[0] = lchi1;
  y[1] = norm_x_2;

  norm_x = aocl_blas_dnrm2( &i_two,
                   y, &i_one );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //          = -sign( chi_1 ) * || x ||_2
  //

  alpha = -dsign( lchi1 ) * norm_x;

  //
  // Overwrite x_2 with u_2:
  //
  //   x_2 := x_2 / ( chi_1 - alpha )
  //

  chi_1_minus_alpha = lchi1 - alpha;

  //
  // Scale X2, chi_1_minus_alpha, alpha
  // if norm factor is very less
  //
  safmin = fla_dlamch("S", 1) / fla_dlamch("E", 1);
  kn = 0;
  if( fabs( chi_1_minus_alpha ) < safmin )
  {
    rsafmn = 1. / safmin;

    for( kn = 1; kn < 20; kn++ )
    {
      aocl_blas_dscal( &m_x2,
               &rsafmn,
               x2, &inc_x2 );
      alpha = alpha * rsafmn;
      chi_1_minus_alpha = chi_1_minus_alpha * rsafmn;
      
      if( fabs( chi_1_minus_alpha ) > safmin )
        break;
    }
  }

  inv_chi_1_minus_alpha = 1.0 / chi_1_minus_alpha;

  aocl_blas_dscal( &m_x2,
          &inv_chi_1_minus_alpha,
          x2, &inc_x2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //        = alpha / ( alpha - chi_1 )
  //

  *tau = -1.0 * alpha * inv_chi_1_minus_alpha;

  //
  // Scale back alpha
  //
  for( ; kn > 0; kn-- )
  {
    alpha = alpha * safmin;
  }

  //
  // Overwrite chi_1 with alpha:
  //
  //   chi_1 := alpha
  //

  *chi_1 = alpha;

  return FLA_SUCCESS;
}

FLA_Error FLA_EXT_dgeqrf( fla_dim_t  m_A, fla_dim_t n_A,
                          double*  buff_A, fla_dim_t cs_A,
                          double*  buff_t,
                          double*  buff_w,
                          fla_dim_t* lwork,
                          fla_dim_t* info )
{
  fla_dim_t min_m_n = fla_min( m_A, n_A );
  fla_dim_t i, rs_A = 1;

  for ( i = 0; i < min_m_n; ++i )
  {
    double* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    double* tau1     = buff_t + i;
    double  alphat   = *alpha11;

    fla_dim_t m_curr   = m_A - i;

    fla_dim_t m_ahead  = m_A - i - 1;
    fla_dim_t n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // Compute Householder transformation for current column
    FLA_EXT_Househ2_l_opd( m_ahead,
                           &alphat,
                           a21, rs_A,
                           tau1 );

    *alpha11 = 1.0;
    if( *tau1 != 0 )
        *tau1 = 1.0 / *tau1;

    // Apply the computed Householder transformation on the matrix
    aocl_lapack_dlarf( "L",
            &m_curr, &n_ahead,
            alpha11, &rs_A,
            tau1,
            a12t, &cs_A,
            buff_w );

    *alpha11 = alphat;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}
