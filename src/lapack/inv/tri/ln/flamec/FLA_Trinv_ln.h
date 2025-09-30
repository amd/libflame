/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Trinv_ln_blk_var1( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_ln_blk_var2( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_ln_blk_var3( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_ln_blk_var4( FLA_Obj A, fla_trinv_t* cntl );

FLA_Error FLA_Trinv_ln_unb_var1( FLA_Obj A );
FLA_Error FLA_Trinv_ln_unb_var2( FLA_Obj A );
FLA_Error FLA_Trinv_ln_unb_var3( FLA_Obj A );
FLA_Error FLA_Trinv_ln_unb_var4( FLA_Obj A );

FLA_Error FLA_Trinv_ln_opt_var1( FLA_Obj A );
FLA_Error FLA_Trinv_ln_ops_var1( fla_dim_t mn_A,
                                 float*    A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opd_var1( fla_dim_t mn_A,
                                 double*   A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opc_var1( fla_dim_t mn_A,
                                 scomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opz_var1( fla_dim_t mn_A,
                                 dcomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );

FLA_Error FLA_Trinv_ln_opt_var2( FLA_Obj A );
FLA_Error FLA_Trinv_ln_ops_var2( fla_dim_t mn_A,
                                 float*    A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opd_var2( fla_dim_t mn_A,
                                 double*   A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opc_var2( fla_dim_t mn_A,
                                 scomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opz_var2( fla_dim_t mn_A,
                                 dcomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );

FLA_Error FLA_Trinv_ln_opt_var3( FLA_Obj A );
FLA_Error FLA_Trinv_ln_ops_var3( fla_dim_t mn_A,
                                 float*    A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opd_var3( fla_dim_t mn_A,
                                 double*   A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opc_var3( fla_dim_t mn_A,
                                 scomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opz_var3( fla_dim_t mn_A,
                                 dcomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );

FLA_Error FLA_Trinv_ln_opt_var4( FLA_Obj A );
FLA_Error FLA_Trinv_ln_ops_var4( fla_dim_t mn_A,
                                 float*    A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opd_var4( fla_dim_t mn_A,
                                 double*   A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opc_var4( fla_dim_t mn_A,
                                 scomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );
FLA_Error FLA_Trinv_ln_opz_var4( fla_dim_t mn_A,
                                 dcomplex* A, fla_dim_t rs_A, fla_dim_t cs_A );
