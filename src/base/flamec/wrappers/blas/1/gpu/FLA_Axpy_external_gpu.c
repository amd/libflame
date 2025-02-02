/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Axpy_external_gpu( FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu )
{
  FLA_Datatype datatype;
  integer          m_B, n_B;
  integer          ldim_A, inc_A;
  integer          ldim_B, inc_B;
  integer          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Axpy_check( alpha, A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  ldim_A   = FLA_Obj_length( A );
  inc_A    = 1;

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  ldim_B   = FLA_Obj_length( B );
  inc_B    = 1;

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_alpha = ( float* ) FLA_FLOAT_PTR( alpha );
    float* buff_A_gpu = ( float* ) A_gpu;
    float* buff_B_gpu = ( float* ) B_gpu;

    for ( i = 0; i < n_B; i++ )
      cublasSaxpy( m_B,
                   *buff_alpha,
                   buff_A_gpu + i * ldim_A, inc_A,
                   buff_B_gpu + i * ldim_B, inc_B );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_alpha = ( double* ) FLA_DOUBLE_PTR( alpha );
    double* buff_A_gpu = ( double* ) A_gpu;
    double* buff_B_gpu = ( double* ) B_gpu;

    for ( i = 0; i < n_B; i++ )
      cublasDaxpy( m_B,
                   *buff_alpha,
                   buff_A_gpu + i * ldim_A, inc_A,
                   buff_B_gpu + i * ldim_B, inc_B );

    break;
  }

  case FLA_COMPLEX:
  {
    cuComplex* buff_alpha = ( cuComplex* ) FLA_COMPLEX_PTR( alpha );
    cuComplex* buff_A_gpu = ( cuComplex* ) A_gpu;
    cuComplex* buff_B_gpu = ( cuComplex* ) B_gpu;

    for ( i = 0; i < n_B; i++ )
      cublasCaxpy( m_B,
                   *buff_alpha,
                   buff_A_gpu + i * ldim_A, inc_A,
                   buff_B_gpu + i * ldim_B, inc_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cuDoubleComplex* buff_alpha = ( cuDoubleComplex* ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    cuDoubleComplex* buff_A_gpu = ( cuDoubleComplex* ) A_gpu;
    cuDoubleComplex* buff_B_gpu = ( cuDoubleComplex* ) B_gpu;

    for ( i = 0; i < n_B; i++ )
      cublasZaxpy( m_B,
                   *buff_alpha,
                   buff_A_gpu + i * ldim_A, inc_A,
                   buff_B_gpu + i * ldim_B, inc_B );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

#endif
