/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_nn_blk_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  FLA_Obj BTL,   BTR,      B00, B01, B02, 
          BBL,   BBR,      B10, B11, B12,
                           B20, B21, B22;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  fla_dim_t b;

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  while ( FLA_Obj_length( BTL ) < FLA_Obj_length( B ) ){

    b = FLA_Determine_blocksize( CR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00, /**/ &B01, &B02,
                        /* ************* */   /* ******************** */
                                                &B10, /**/ &B11, &B12,
                           BBL, /**/ BBR,       &B20, /**/ &B21, &B22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );

    // Loop Invariant:
    // CL = 
    // CR = 

    /*------------------------------------------------------------*/

    // C1 = sylv( A, B11, C1 );
    FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       isgn, A, B11, C1, scale,
                       FLA_Cntl_sub_sylv1( cntl ) );

    // C2 = C2 -/+ C1 * B12;
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_NEGATE( isgn ), C1, B12, FLA_ONE, C2,
                       FLA_Cntl_sub_gemm1( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00, B01, /**/ B02,
                                                     B10, B11, /**/ B12,
                            /* ************** */  /* ****************** */
                              &BBL, /**/ &BBR,       B20, B21, /**/ B22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}
