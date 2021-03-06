/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Tevd_francis_v_opt_var1() -------------------------------------------

FLA_Error FLA_Tevd_francis_v_opt_var1( FLA_Obj shift, FLA_Obj g, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tevd_francis_v_ops_var1( integer       m_A,
                                       float*    buff_shift,
                                       scomplex* buff_g, integer inc_g, 
                                       float*    buff_d, integer inc_d, 
                                       float*    buff_e, integer inc_e ); 
FLA_Error FLA_Tevd_francis_v_opd_var1( integer       m_A,
                                       double*   buff_shift,
                                       dcomplex* buff_g, integer inc_g, 
                                       double*   buff_d, integer inc_d, 
                                       double*   buff_e, integer inc_e ); 

