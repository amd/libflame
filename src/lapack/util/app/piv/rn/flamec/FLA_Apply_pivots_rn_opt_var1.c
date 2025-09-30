/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_pivots_rn_opt_var1( FLA_Obj p, FLA_Obj A )
{
  FLA_Datatype datatype;
  fla_dim_t          m_A;
  fla_dim_t          rs_A, cs_A;
  fla_dim_t          inc_p;
  fla_dim_t          k1_0, k2_0;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );

  // Swap the stride; FLA_Apply_pivots_ln_ops_var1 already consider the memory access pattern.
  cs_A     = FLA_Obj_row_stride( A );
  rs_A     = FLA_Obj_col_stride( A );

  inc_p    = -FLA_Obj_vector_inc( p );

  // Use zero-based indices.
  k1_0     = 0;
  k2_0     = ( fla_dim_t ) FLA_Obj_vector_dim( p ) - 1;

  switch ( datatype )
  {
    case FLA_INT:
    {
      fla_dim_t*   buff_A = FLA_INT_PTR( A );
      fla_dim_t*   buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opi_var1( m_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      fla_dim_t*   buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_ops_var1( m_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      fla_dim_t*    buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opd_var1( m_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      fla_dim_t*      buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opc_var1( m_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      fla_dim_t*      buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opz_var1( m_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }
  }

  return FLA_SUCCESS;
}

