/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Tevd_eigval_n_opt_var1() --------------------------------------------

FLA_Error FLA_Tevd_eigval_n_opt_var1( FLA_Obj G, FLA_Obj d, FLA_Obj e, FLA_Obj n_iter );
FLA_Error FLA_Tevd_eigval_n_ops_var1( fla_dim_t       m_A,
                                      fla_dim_t       n_G,
                                      float*    buff_d, fla_dim_t inc_d, 
                                      float*    buff_e, fla_dim_t inc_e,
                                      fla_dim_t*      n_iter );
FLA_Error FLA_Tevd_eigval_n_opd_var1( fla_dim_t       m_A,
                                      fla_dim_t       n_G,
                                      double*   buff_d, fla_dim_t inc_d, 
                                      double*   buff_e, fla_dim_t inc_e,
                                      fla_dim_t*      n_iter );

