/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Utility-level BLAS-like prototypes --------------------------------------

// --- constant-generating functions ---

float    bl1_s2( void );
double   bl1_d2( void );
scomplex bl1_c2( void );
dcomplex bl1_z2( void );
float    bl1_s1( void );
double   bl1_d1( void );
scomplex bl1_c1( void );
dcomplex bl1_z1( void );
float    bl1_s1h( void );
double   bl1_d1h( void );
scomplex bl1_c1h( void );
dcomplex bl1_z1h( void );
float    bl1_s0( void );
double   bl1_d0( void );
scomplex bl1_c0( void );
dcomplex bl1_z0( void );
float    bl1_sm1h( void );
double   bl1_dm1h( void );
scomplex bl1_cm1h( void );
dcomplex bl1_zm1h( void );
float    bl1_sm1( void );
double   bl1_dm1( void );
scomplex bl1_cm1( void );
dcomplex bl1_zm1( void );
float    bl1_sm2( void );
double   bl1_dm2( void );
scomplex bl1_cm2( void );
dcomplex bl1_zm2( void );

// --- allocv ---

void*     bl1_vallocv( uinteger n_elem, uinteger elem_size );
fla_dim_t*      bl1_iallocv( uinteger n_elem );
float*    bl1_sallocv( uinteger n_elem );
double*   bl1_dallocv( uinteger n_elem );
scomplex* bl1_callocv( uinteger n_elem );
dcomplex* bl1_zallocv( uinteger n_elem );

// --- allocm ---

void*     bl1_vallocm( uinteger m, uinteger n, uinteger elem_size );
fla_dim_t*      bl1_iallocm( uinteger m, uinteger n );
float*    bl1_sallocm( uinteger m, uinteger n );
double*   bl1_dallocm( uinteger m, uinteger n );
scomplex* bl1_callocm( uinteger m, uinteger n );
dcomplex* bl1_zallocm( uinteger m, uinteger n );

// --- apdiagmv ---

void bl1_sapdiagmv( side1_t side, conj1_t conj, fla_dim_t m, fla_dim_t n, float*    x, fla_dim_t incx, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dapdiagmv( side1_t side, conj1_t conj, fla_dim_t m, fla_dim_t n, double*   x, fla_dim_t incx, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csapdiagmv( side1_t side, conj1_t conj, fla_dim_t m, fla_dim_t n, float*    x, fla_dim_t incx, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_capdiagmv( side1_t side, conj1_t conj, fla_dim_t m, fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zdapdiagmv( side1_t side, conj1_t conj, fla_dim_t m, fla_dim_t n, double*   x, fla_dim_t incx, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zapdiagmv( side1_t side, conj1_t conj, fla_dim_t m, fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- create_contigm ---

void bl1_screate_contigm( fla_dim_t m, fla_dim_t n, float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dcreate_contigm( fla_dim_t m, fla_dim_t n, double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_ccreate_contigm( fla_dim_t m, fla_dim_t n, scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zcreate_contigm( fla_dim_t m, fla_dim_t n, dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- create_contigmt ---

void bl1_screate_contigmt( trans1_t trans_dims, fla_dim_t m, fla_dim_t n, float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dcreate_contigmt( trans1_t trans_dims, fla_dim_t m, fla_dim_t n, double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_ccreate_contigmt( trans1_t trans_dims, fla_dim_t m, fla_dim_t n, scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zcreate_contigmt( trans1_t trans_dims, fla_dim_t m, fla_dim_t n, dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- create_contigmr ---

void bl1_screate_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dcreate_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_ccreate_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zcreate_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- create_contigmsr ---

void bl1_screate_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dcreate_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_ccreate_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zcreate_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- free_contigm ---

void bl1_sfree_contigm( float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dfree_contigm( double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_cfree_contigm( scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zfree_contigm( dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- free_saved_contigm ---

void bl1_sfree_saved_contigm( fla_dim_t m, fla_dim_t n, float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dfree_saved_contigm( fla_dim_t m, fla_dim_t n, double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_cfree_saved_contigm( fla_dim_t m, fla_dim_t n, scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zfree_saved_contigm( fla_dim_t m, fla_dim_t n, dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- free_saved_contigmr ---

void bl1_sfree_saved_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dfree_saved_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_cfree_saved_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zfree_saved_contigmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- free_saved_contigmsr ---

void bl1_sfree_saved_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, float**    a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_dfree_saved_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, double**   a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_cfree_saved_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, scomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );
void bl1_zfree_saved_contigmsr( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a_save, fla_dim_t a_rs_save, fla_dim_t a_cs_save, dcomplex** a, fla_dim_t* a_rs, fla_dim_t* a_cs );

// --- ewinvscalv ---

void bl1_sewinvscalv( conj1_t conj, fla_dim_t n, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_dewinvscalv( conj1_t conj, fla_dim_t n, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_csewinvscalv( conj1_t conj, fla_dim_t n, float*    x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_cewinvscalv( conj1_t conj, fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zdewinvscalv( conj1_t conj, fla_dim_t n, double*   x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );
void bl1_zewinvscalv( conj1_t conj, fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

// --- ewscalmt ---

void bl1_sewinvscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dewinvscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_csewinvscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cewinvscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zdewinvscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zewinvscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- ewscalv ---

void bl1_sewscalv( conj1_t conj, fla_dim_t n, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_dewscalv( conj1_t conj, fla_dim_t n, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_csewscalv( conj1_t conj, fla_dim_t n, float*    x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_cewscalv( conj1_t conj, fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zdewscalv( conj1_t conj, fla_dim_t n, double*   x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );
void bl1_zewscalv( conj1_t conj, fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

// --- ewscalmt ---

void bl1_sewscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dewscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_csewscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cewscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zdewscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zewscalmt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- free ---

void bl1_vfree( void*     p );
void bl1_ifree( fla_dim_t*      p );
void bl1_sfree( float*    p );
void bl1_dfree( double*   p );
void bl1_cfree( scomplex* p );
void bl1_zfree( dcomplex* p );

// --- inverts ---

void bl1_sinverts( conj1_t conj, float*    alpha );
void bl1_dinverts( conj1_t conj, double*   alpha );
void bl1_cinverts( conj1_t conj, scomplex* alpha );
void bl1_zinverts( conj1_t conj, dcomplex* alpha );

// --- invert2s ---

void bl1_sinvert2s( conj1_t conj, float*    alpha, float*    beta );
void bl1_dinvert2s( conj1_t conj, double*   alpha, double*   beta );
void bl1_cinvert2s( conj1_t conj, scomplex* alpha, scomplex* beta );
void bl1_zinvert2s( conj1_t conj, dcomplex* alpha, dcomplex* beta );

// --- invertv ---

void bl1_sinvertv( conj1_t conj, fla_dim_t n, float*    x, fla_dim_t incx );
void bl1_dinvertv( conj1_t conj, fla_dim_t n, double*   x, fla_dim_t incx );
void bl1_cinvertv( conj1_t conj, fla_dim_t n, scomplex* x, fla_dim_t incx );
void bl1_zinvertv( conj1_t conj, fla_dim_t n, dcomplex* x, fla_dim_t incx );

// --- ident ---

void bl1_sident( fla_dim_t m, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dident( fla_dim_t m, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cident( fla_dim_t m, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zident( fla_dim_t m, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- maxabsv ---

void bl1_smaxabsv( fla_dim_t n, float*    x, fla_dim_t incx, float*  maxabs );
void bl1_dmaxabsv( fla_dim_t n, double*   x, fla_dim_t incx, double* maxabs );
void bl1_cmaxabsv( fla_dim_t n, scomplex* x, fla_dim_t incx, float*  maxabs );
void bl1_zmaxabsv( fla_dim_t n, dcomplex* x, fla_dim_t incx, double* maxabs );

// --- maxabsm ---

void bl1_smaxabsm( fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*  maxabs );
void bl1_dmaxabsm( fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double* maxabs );
void bl1_cmaxabsm( fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*  maxabs );
void bl1_zmaxabsm( fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double* maxabs );

// --- maxabsmr ---

void bl1_smaxabsmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*  maxabs );
void bl1_dmaxabsmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double* maxabs );
void bl1_cmaxabsmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*  maxabs );
void bl1_zmaxabsmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double* maxabs );

// --- rands ---

void bl1_srands( float*    alpha );
void bl1_drands( double*   alpha );
void bl1_crands( scomplex* alpha );
void bl1_zrands( dcomplex* alpha );

// --- randv ---

void bl1_srandv( fla_dim_t n, float*    x, fla_dim_t incx );
void bl1_drandv( fla_dim_t n, double*   x, fla_dim_t incx );
void bl1_crandv( fla_dim_t n, scomplex* x, fla_dim_t incx );
void bl1_zrandv( fla_dim_t n, dcomplex* x, fla_dim_t incx );

// --- randm ---

void bl1_srandm( fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_drandm( fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_crandm( fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zrandm( fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- randmr ---
void bl1_srandmr( uplo1_t uplo, diag1_t diag, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_drandmr( uplo1_t uplo, diag1_t diag, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_crandmr( uplo1_t uplo, diag1_t diag, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zrandmr( uplo1_t uplo, diag1_t diag, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- set_contig_strides ---

void bl1_set_contig_strides( fla_dim_t m, fla_dim_t n, fla_dim_t* rs, fla_dim_t* cs );

// --- set_dims_with_side ---

void bl1_set_dim_with_side( side1_t side, fla_dim_t m, fla_dim_t n, fla_dim_t* dim_new );

// --- set_dims_with_trans ---

void bl1_set_dims_with_trans( trans1_t trans, fla_dim_t m, fla_dim_t n, fla_dim_t* m_new, fla_dim_t* n_new );

// --- setv ---

void bl1_isetv( fla_dim_t m, fla_dim_t*      sigma, fla_dim_t*      x, fla_dim_t incx );
void bl1_ssetv( fla_dim_t m, float*    sigma, float*    x, fla_dim_t incx );
void bl1_dsetv( fla_dim_t m, double*   sigma, double*   x, fla_dim_t incx );
void bl1_csetv( fla_dim_t m, scomplex* sigma, scomplex* x, fla_dim_t incx );
void bl1_zsetv( fla_dim_t m, dcomplex* sigma, dcomplex* x, fla_dim_t incx );

// --- setm ---

void bl1_isetm( fla_dim_t m, fla_dim_t n, fla_dim_t*      sigma, fla_dim_t*      a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_ssetm( fla_dim_t m, fla_dim_t n, float*    sigma, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dsetm( fla_dim_t m, fla_dim_t n, double*   sigma, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csetm( fla_dim_t m, fla_dim_t n, scomplex* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zsetm( fla_dim_t m, fla_dim_t n, dcomplex* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- setmr ---

void bl1_ssetmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    sigma, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dsetmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   sigma, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csetmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zsetmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- setdiag ---

void bl1_isetdiag( fla_dim_t offset, fla_dim_t m, fla_dim_t n, fla_dim_t*      sigma, fla_dim_t*      a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_ssetdiag( fla_dim_t offset, fla_dim_t m, fla_dim_t n, float*    sigma, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dsetdiag( fla_dim_t offset, fla_dim_t m, fla_dim_t n, double*   sigma, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csetdiag( fla_dim_t offset, fla_dim_t m, fla_dim_t n, scomplex* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zsetdiag( fla_dim_t offset, fla_dim_t m, fla_dim_t n, dcomplex* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- scalediag ---

void bl1_sscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, float*    sigma, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, double*   sigma, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, scomplex* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, dcomplex* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, float*    sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zdscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, double*   sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- shiftdiag ---

void bl1_sshiftdiag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, float*    sigma, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dshiftdiag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, double*   sigma, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cshiftdiag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, scomplex* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zshiftdiag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, dcomplex* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csshiftdiag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, float*    sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zdshiftdiag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, double*   sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- symmize ---

void bl1_ssymmize( conj1_t conj, uplo1_t uplo, fla_dim_t m, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dsymmize( conj1_t conj, uplo1_t uplo, fla_dim_t m, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csymmize( conj1_t conj, uplo1_t uplo, fla_dim_t m, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zsymmize( conj1_t conj, uplo1_t uplo, fla_dim_t m, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

