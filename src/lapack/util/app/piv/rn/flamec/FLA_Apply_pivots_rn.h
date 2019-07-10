/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_pivots_rn_opt_var1( FLA_Obj p, FLA_Obj A );
FLA_Error FLA_Apply_pivots_rn_ops_var1( fla_dim_t n, 
                                        float*    a, fla_dim_t a_rs, fla_dim_t a_cs, 
                                        fla_dim_t k1, 
                                        fla_dim_t k2, 
                                        fla_dim_t* p, fla_dim_t incp );
FLA_Error FLA_Apply_pivots_rn_opd_var1( fla_dim_t n, 
                                        double*   a, fla_dim_t a_rs, fla_dim_t a_cs, 
                                        fla_dim_t k1, 
                                        fla_dim_t k2, 
                                        fla_dim_t* p, fla_dim_t incp );
FLA_Error FLA_Apply_pivots_rn_opc_var1( fla_dim_t n, 
                                        scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, 
                                        fla_dim_t k1, 
                                        fla_dim_t k2, 
                                        fla_dim_t* p, fla_dim_t incp );
FLA_Error FLA_Apply_pivots_rn_opz_var1( fla_dim_t n, 
                                        dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, 
                                        fla_dim_t k1, 
                                        fla_dim_t k2, 
                                        fla_dim_t* p, fla_dim_t incp );
