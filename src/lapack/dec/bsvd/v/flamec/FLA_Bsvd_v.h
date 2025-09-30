/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Bsvd_iteracc_v.h"
#include "FLA_Bsvd_sinval_v.h"
#include "FLA_Bsvd_francis_v.h"

// --- FLA_Bsvd_compute_shift() ------------------------------------------------

FLA_Error FLA_Bsvd_compute_shift( FLA_Obj tol, FLA_Obj sminl, FLA_Obj smax, FLA_Obj d, FLA_Obj e, FLA_Obj shift );
FLA_Error FLA_Bsvd_compute_shift_ops( fla_dim_t       m_A,
                                      float     tol,
                                      float     sminl,
                                      float     smax,
                                      float*    buff_d, fla_dim_t inc_d,
                                      float*    buff_e, fla_dim_t inc_e,
                                      float*    shift );
FLA_Error FLA_Bsvd_compute_shift_opd( fla_dim_t       m_A,
                                      double    tol,
                                      double    sminl,
                                      double    smax,
                                      double*   buff_d, fla_dim_t inc_d,
                                      double*   buff_e, fla_dim_t inc_e,
                                      double*   shift );

// --- FLA_Bsvd_compute_tol_thresh() -------------------------------------------

FLA_Error FLA_Bsvd_compute_tol_thresh( FLA_Obj tolmul, FLA_Obj maxit, FLA_Obj d, FLA_Obj e, FLA_Obj tol, FLA_Obj thresh );
FLA_Error FLA_Bsvd_compute_tol_thresh_ops( fla_dim_t       m_A,
                                           float     tolmul,
                                           float     maxit,
                                           float*    buff_d, fla_dim_t inc_d, 
                                           float*    buff_e, fla_dim_t inc_e, 
                                           float*    tol,
                                           float*    thresh );
FLA_Error FLA_Bsvd_compute_tol_thresh_opd( fla_dim_t       m_A,
                                           double    tolmul,
                                           double    maxit,
                                           double*   buff_d, fla_dim_t inc_d, 
                                           double*   buff_e, fla_dim_t inc_e, 
                                           double*   tol,
                                           double*   thresh );

// --- FLA_Bsvd_find_converged() -----------------------------------------------

FLA_Error FLA_Bsvd_find_converged( FLA_Obj tol, FLA_Obj d, FLA_Obj e, FLA_Obj sminl );
FLA_Error FLA_Bsvd_find_converged_ops( fla_dim_t       m_A,
                                       float     tol, 
                                       float*    buff_d, fla_dim_t inc_d, 
                                       float*    buff_e, fla_dim_t inc_e,
                                       float*    sminl );
FLA_Error FLA_Bsvd_find_converged_opd( fla_dim_t       m_A,
                                       double    tol, 
                                       double*   buff_d, fla_dim_t inc_d, 
                                       double*   buff_e, fla_dim_t inc_e,
                                       double*   sminl );

// --- FLA_Bsvd_find_max_fla_min() -------------------------------------------------

FLA_Error FLA_Bsvd_find_max_fla_min( FLA_Obj d, FLA_Obj e, FLA_Obj smax, FLA_Obj smin );
FLA_Error FLA_Bsvd_find_max_min_ops( fla_dim_t       m_A,
                                     float*    buff_d, fla_dim_t inc_d, 
                                     float*    buff_e, fla_dim_t inc_e, 
                                     float*    smax,
                                     float*    smin );
FLA_Error FLA_Bsvd_find_max_min_opd( fla_dim_t       m_A,
                                     double*   buff_d, fla_dim_t inc_d, 
                                     double*   buff_e, fla_dim_t inc_e, 
                                     double*   smax,
                                     double*   smin );

// --- FLA_Bsvd_find_submatrix() -----------------------------------------------

FLA_Error FLA_Bsvd_find_submatrix_ops( fla_dim_t       mn_A,
                                       fla_dim_t       ij_begin,
                                       float*    buff_d, fla_dim_t inc_d,
                                       float*    buff_e, fla_dim_t inc_e,
                                       fla_dim_t*      ijTL,
                                       fla_dim_t*      ijBR );
FLA_Error FLA_Bsvd_find_submatrix_opd( fla_dim_t       mn_A,
                                       fla_dim_t       ij_begin,
                                       double*   buff_d, fla_dim_t inc_d,
                                       double*   buff_e, fla_dim_t inc_e,
                                       fla_dim_t*      ijTL,
                                       fla_dim_t*      ijBR );

// --- FLA_Bsvd_v_opt_var1() ---------------------------------------------------

FLA_Error FLA_Bsvd_v_opt_var1( fla_dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, FLA_Obj U, FLA_Obj V, fla_dim_t b_alg );
FLA_Error FLA_Bsvd_v_ops_var1( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               scomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               float*    buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               float*    buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_v_opd_var1( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               dcomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               double*   buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               double*   buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_v_opc_var1( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               scomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               scomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               scomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_v_opz_var1( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               dcomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               dcomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               dcomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );

// --- FLA_Bsvd_v_opt_var2() ---------------------------------------------------

FLA_Error FLA_Bsvd_v_opt_var2( fla_dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, FLA_Obj RG, FLA_Obj RH, FLA_Obj W, FLA_Obj U, FLA_Obj V, fla_dim_t b_alg );
FLA_Error FLA_Bsvd_v_ops_var2( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               scomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               float*    buff_RG, fla_dim_t rs_RG, fla_dim_t cs_RG,
                               float*    buff_RH, fla_dim_t rs_RH, fla_dim_t cs_RH,
                               float*    buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               float*    buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               float*    buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_v_opd_var2( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               dcomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               double*   buff_RG, fla_dim_t rs_RG, fla_dim_t cs_RG,
                               double*   buff_RH, fla_dim_t rs_RH, fla_dim_t cs_RH,
                               double*   buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               double*   buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               double*   buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_v_opc_var2( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               float*    buff_d, fla_dim_t inc_d, 
                               float*    buff_e, fla_dim_t inc_e,
                               scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               scomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               float*    buff_RG, fla_dim_t rs_RG, fla_dim_t cs_RG,
                               float*    buff_RH, fla_dim_t rs_RH, fla_dim_t cs_RH,
                               scomplex* buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               scomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               scomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_v_opz_var2( fla_dim_t       min_m_n,
                               fla_dim_t       m_U,
                               fla_dim_t       m_V,
                               fla_dim_t       n_GH,
                               fla_dim_t       n_iter_max,
                               double*   buff_d, fla_dim_t inc_d, 
                               double*   buff_e, fla_dim_t inc_e,
                               dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                               dcomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                               double*   buff_RG, fla_dim_t rs_RG, fla_dim_t cs_RG,
                               double*   buff_RH, fla_dim_t rs_RH, fla_dim_t cs_RH,
                               dcomplex* buff_W, fla_dim_t rs_W, fla_dim_t cs_W,
                               dcomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                               dcomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                               fla_dim_t       b_alg );

