/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
/*
 *  Modifications Copyright (c) 2024 Advanced Micro Devices, Inc.  All rights reserved.
 */

FLA_Error FLA_Svd_ext_u_unb_var1( FLA_Svd_type jobu, FLA_Svd_type jobv, 
                                  fla_dim_t n_iter_max,
                                  FLA_Obj A, FLA_Obj s, FLA_Obj V, FLA_Obj U,
                                  fla_dim_t k_accum,
                                  fla_dim_t b_alg );
int lapack_dbdsqr(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *
	          nru, aocl_int64_t *ncc, doublereal *d__, doublereal *e, doublereal *vt, 
	          aocl_int64_t *ldvt, doublereal *u, aocl_int64_t *ldu, doublereal *c__, aocl_int64_t *
	          ldc, doublereal *work, aocl_int64_t *info);
int lapack_dbdsqr_small(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru,
                                doublereal *d__, doublereal *e,
                                doublereal *vt, aocl_int64_t *ldvt,
                                doublereal *u, aocl_int64_t *ldu,
                                aocl_int64_t *info);
int lapack_dgebd2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *
	          lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	          taup, doublereal *work, aocl_int64_t *info);
int lapack_dgebrd(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *
	          lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	          taup, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dgelqf(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *
	          lda, doublereal *tau, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dgelq2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t * lda, 
                  doublereal *tau, doublereal *work, aocl_int64_t *info);
int lapack_dgesvd(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, 
	          doublereal *a, aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *
	          ldu, doublereal *vt, aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *lwork, 
	          aocl_int64_t *info);
int lapack_dorg2r(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *
	          a, aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *info);
int lapack_dorgbr(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, 
	          doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *work, 
	          aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dorgl2(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *
	          a, aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *info);
int lapack_dorglq(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *
	          a, aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *lwork, 
	          aocl_int64_t *info);
int lapack_dorgqr(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *
	          a, aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *lwork, 
	          aocl_int64_t *info);
int lapack_dorm2r(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, 
                  doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal * c__, 
                  aocl_int64_t *ldc, doublereal *work, aocl_int64_t *info);
int lapack_dormbr(char *vect, char *side, char *trans, aocl_int64_t *m, 
	          aocl_int64_t *n, aocl_int64_t *k, doublereal *a, aocl_int64_t *lda, doublereal *tau, 
	          doublereal *c__, aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, 
	          aocl_int64_t *info);
int lapack_dormlq(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, 
                  aocl_int64_t *k, doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *
	          c__, aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dorml2(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, 
                  doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal * c__, 
                  aocl_int64_t *ldc, doublereal *work, aocl_int64_t *info); 
int lapack_dormqr(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, 
	          aocl_int64_t *k, doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *
	          c__, aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int  dgesvd2x2(   char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n,
                  doublereal *a, aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *
                  ldu, doublereal *vt, aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);



