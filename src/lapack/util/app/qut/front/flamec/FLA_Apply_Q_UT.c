/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_apqut_t* fla_apqut_cntl_leaf;
extern TLS_CLASS_SPEC fla_apqut_t* fla_apqut_cntl;

FLA_Error FLA_Apply_Q_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_Q_UT_check( side, trans, direct, storev, A, T, W, B );

  // Invoke FLA_Apply_Q_UT_internal() with the standard control tree.
  r_val = FLA_Apply_Q_UT_internal( side, trans, direct, storev, A, T, W, B, fla_apqut_cntl_leaf );

  return r_val;
}

