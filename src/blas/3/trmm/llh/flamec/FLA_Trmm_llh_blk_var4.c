/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trmm_llh_blk_var4( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  fla_dim_t b;

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  while ( FLA_Obj_width( BR ) < FLA_Obj_width( B ) ){

    b = FLA_Determine_blocksize( BL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &B1, /**/ &B2,
                           b, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* B1 = tril( A )' *  B1 */
    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diagA,
                       alpha, A, B1,
                       FLA_Cntl_sub_trmm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ B1, B2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
