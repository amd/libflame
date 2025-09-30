/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Tevd_iteracc_v.h"
#include "FLA_Tevd_eigval_v.h"
#include "FLA_Tevd_francis_v.h"

// --- FLA_Tevd_compute_scaling() ----------------------------------------------

FLA_Error FLA_Tevd_compute_scaling_ops( fla_dim_t       m_A,
                                        float*    buff_d, fla_dim_t inc_d, 
                                        float*    buff_e, fla_dim_t inc_e,
                                        float*    sigma );
FLA_Error FLA_Tevd_compute_scaling_opd( fla_dim_t       m_A,
                                        double*   buff_d, fla_dim_t inc_d, 
                                        double*   buff_e, fla_dim_t inc_e,
                                        double*   sigma );

// --- FLA_Tevd_find_submatrix() -----------------------------------------------

FLA_Error FLA_Tevd_find_submatrix_ops( fla_dim_t       m_A,
                                       fla_dim_t       ij_begin,
                                       float*    buff_d, fla_dim_t inc_d, 
                                       float*    buff_e, fla_dim_t inc_e,
                                       fla_dim_t*      ijTL,
                                       fla_dim_t*      ijBR );
FLA_Error FLA_Tevd_find_submatrix_opd( fla_dim_t       m_A,
                                       fla_dim_t       ij_begin,
                                       double*   buff_d, fla_dim_t inc_d, 
                                       double*   buff_e, fla_dim_t inc_e,
                                       fla_dim_t*      ijTL,
                                       fla_dim_t*      ijBR );

// --- FLA_Tevd_find_perfshift() -----------------------------------------------

FLA_Error FLA_Tevd_find_perfshift_ops( fla_dim_t       m_d,
                                       fla_dim_t       m_l,
                                       float*    buff_d, fla_dim_t inc_d, 
                                       float*    buff_e, fla_dim_t inc_e, 
                                       float*    buff_l, fla_dim_t inc_l, 
                                       fla_dim_t*      buff_lstat, fla_dim_t inc_lstat, 
                                       float*    buff_pu, fla_dim_t inc_pu, 
                                       fla_dim_t*      ij_shift );
FLA_Error FLA_Tevd_find_perfshift_opd( fla_dim_t       m_d,
                                       fla_dim_t       m_l,
                                       double*   buff_d, fla_dim_t inc_d, 
                                       double*   buff_e, fla_dim_t inc_e, 
                                       double*   buff_l, fla_dim_t inc_l, 
                                       fla_dim_t*      buff_lstat, fla_dim_t inc_lstat, 
                                       double*   buff_pu, fla_dim_t inc_pu, 
                                       fla_dim_t*      ij_shift );

// --- FLA_Norm1_tridiag() -----------------------------------------------------

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm );
FLA_Error FLA_Norm1_tridiag_ops( fla_dim_t       m_A,
                                 float*    buff_d, fla_dim_t inc_d, 
                                 float*    buff_e, fla_dim_t inc_e,
                                 float*    norm );
FLA_Error FLA_Norm1_tridiag_opd( fla_dim_t       m_A,
                                 double*   buff_d, fla_dim_t inc_d, 
                                 double*   buff_e, fla_dim_t inc_e,
                                 double*   norm );

// --- FLA_Tevd_v_opt_var1() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var1( fla_dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U, fla_dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var1( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_iter_max,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               float*    buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );
FLA_Error FLA_Tevd_v_opd_var1( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_iter_max,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               double*   buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );
FLA_Error FLA_Tevd_v_opc_var1( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_iter_max,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               scomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );
FLA_Error FLA_Tevd_v_opz_var1( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_iter_max,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               dcomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );

// --- FLA_Tevd_v_opt_var2() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var2( fla_dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj U, fla_dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var2( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_G_extra,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               float*    buff_R, fla_dim_t rs_R, fla_dim_t cs_R,
                               float*    buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               float*    buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );
FLA_Error FLA_Tevd_v_opd_var2( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_G_extra,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               double*   buff_R, fla_dim_t rs_R, fla_dim_t cs_R,
                               double*   buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               double*   buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );
FLA_Error FLA_Tevd_v_opc_var2( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_G_extra,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               float*    buff_R, fla_dim_t rs_R, fla_dim_t cs_R,
                               scomplex* buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               scomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );
FLA_Error FLA_Tevd_v_opz_var2( fla_dim_t       m_A,
                               fla_dim_t       m_U,
                               fla_dim_t       n_G,
                               fla_dim_t       n_G_extra,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               double*   buff_R, fla_dim_t rs_R, fla_dim_t cs_R,
                               dcomplex* buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               dcomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               fla_dim_t       b_alg );

