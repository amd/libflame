/******************************************************************************
 * Copyright (C) 2021-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface.hh
 *  libflame_interface.hh defines all the Libflame CPP templated public
 *  interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_HH
#define LIBFLAME_INTERFACE_HH

#include "libflame.hh"

	/** @defgroup AOCL_LAPACK AOCL-LAPACK
	* @{
	*/

	/** @defgroup LAPACK LAPACK
	* @ingroup AOCL_LAPACK
	* @{
	*/

#include "libflame_interface_lin.hh"
#include "libflame_interface_lsq.hh"
#include "libflame_interface_svd.hh"
#include "libflame_interface_eig.hh"
#include "libflame_interface_orth.hh"
#include "libflame_interface_aux.hh"
#include "libflame_interface_blas.hh"

	/** @}*/ // end of LAPACK
	/** @} */ // end of AOCL_LAPACK

namespace libflame
{  
   // Missing API interfaces added without doxygen documentation

	 template<typename T>
	 void geqpf(integer* m, integer* n, T* a, integer* lda, 
               integer* jpvt, T* tau, T* work, integer* info)
	 {
	   geqpf(m, n, a, lda, jpvt, tau, work, info);
	 }
	 
	 template<typename T, typename Ta>
	 void geqpf(integer* m, integer* n, T* a, integer* lda, integer* jpvt, 
               T* tau, T* work, Ta* rwork, integer* info)
	 {
	   geqpf(m, n, a, lda, jpvt, tau, work, rwork, info);
	 }

	 template< typename T >
	 void syr(char* uplo, integer* n, T alpha,  T* x, integer* incx, T* a, integer* lda)
	 {
	   syr(uplo, n, alpha, x, incx, a, lda);
	 }

	 inline void lag2d(integer* m, integer* n, float* sa, integer* ldsa, 
                      double* a, integer* lda, integer* info)
	 {
	   slag2d(m, n, sa, ldsa, a, lda, info);
	 }

	 inline void lag2s(integer*m, integer* n, double* a, integer* lda, float* sa,
                      integer *ldsa, integer* info)
	 {
	   dlag2s(m, n, a, lda, sa, ldsa, info);
	 }

	 inline void lag2z(integer*m, integer* n, scomplex* sa, integer*ldsa, dcomplex* a,
                      integer* lda, integer* info)
	 {
	   clag2z(m, n, sa, ldsa, a, lda, info);
	 }

	 inline void lag2c(integer*m, integer* n, dcomplex* a, integer*lda, scomplex* sa, 
                      integer* ldsa, integer* info)
	 {
	   zlag2c(m, n, a, lda, sa, ldsa, info);
	 }

	 template< typename T >
	 void gelsx(integer* m, integer* n, integer* nrhs, T* a, integer* lda, T* b, 
               integer* ldb, integer* jpvt, T rcond, integer* rank, T* work, integer* info)
	 {
	   gelsx(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, info);
	 }
	 template< typename T , typename Ta >
	 void gelsx(integer* m, integer* n, integer* nrhs, T* a, integer* lda, T* b, 
               integer* ldb, integer* jpvt, Ta rcond, integer* rank, T* work, Ta* rwork, integer* info)
	 {
	   gelsx(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, rwork, info);
	 }

	 inline void sgesv(integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv,
                      double* b, integer* ldb, double* x, integer* ldx, double* work, float* swork, 
                      integer* iter, integer* info)
	 {
	   dsgesv_(n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, iter, info);
	 }

	 inline void cgesv(integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b,
                      integer* ldb, dcomplex* x, integer* ldx, dcomplex *work, scomplex *swork, double *rwork,
                      integer *iter, integer *info)
	 {
	   zcgesv_(n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, rwork, iter, info);
	 }

	 template< typename T >
	 void ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k,
               integer* l, T* a, integer* lda, T* b, integer* ldb, T* alpha, T* beta, T* u, integer* ldu,
               T* v, integer* ldv, T* q, integer* ldq, T* work, integer* iwork, integer* info)
	 {
	   ggsvd(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work,
            iwork, info);
	 }
	 template< typename T , typename Ta >
	 void ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, 
               integer* l, T* a, integer* lda, T* b, integer* ldb, Ta* alpha, Ta* beta, T* u, integer* ldu, 
               T* v, integer* ldv, T* q, integer* ldq, T* work, Ta* rwork, integer* iwork, integer* info)
	 {
	   ggsvd(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work,
      rwork, iwork, info);
	 }

	 template< typename T >
	 void ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, T* a, integer* lda,
               T* b, integer* ldb, T* tola, T* tolb, integer* k, integer* l, T* u, integer* ldu, T* v, 
               integer* ldv, T* q, integer* ldq, integer* iwork, T* tau, T* work, integer* info)
	 {
	   ggsvp(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork,
            tau, work, info);
	 }
	 template< typename T, typename Ta >
	 void ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, T* a, integer* lda,
               T* b, integer* ldb, Ta* tola, Ta* tolb, integer* k, integer* l, T* u, integer* ldu, T* v, 
               integer* ldv, T* q, integer* ldq, integer* iwork, Ta* rwork, T* tau, T* work, integer* info)
	 {
	   ggsvp(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, 
            rwork, tau, work, info);
    }

    template< typename T >
	 void lartgs(T* x, T* y, T* sigma, T* cs, T* sn)
	 {
	   lartgs(x, y, sigma, cs, sn);
	 }

	 template< typename T >
	 void pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  T* ab, integer* ldab,  T* afb,
               integer* ldafb,  T* b, integer* ldb, T* x, integer* ldx, T* ferr, T* berr, T* work, 
               integer* iwork, integer* info)
	 {
	   pbrfs(uplo, n, kd, nrhs,  ab, ldab,  afb, ldafb,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
	 }
	 template< typename T, typename Ta >
	 void pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  T* ab, integer* ldab,  T* afb, 
               integer* ldafb,  T* b, integer* ldb, T* x, integer* ldx, Ta* ferr, Ta* berr, T* work, 
               Ta* rwork, integer* info)
	 {
	   pbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, rwork, info);
	 }

	 template< typename T >
	 void pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, T* ab, integer* ldab, T* b, integer* ldb, 
              integer* info)
	 {
	   pbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
	 }

	 template< typename T >
	 void pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, T* ab, integer* ldab, 
               T* afb, integer* ldafb, char* equed, T* s, T* b, integer* ldb, T* x, integer* ldx, T* rcond, 
               T* ferr, T* berr, T* work, integer* iwork, integer* info)
	 {
	   pbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work,
            iwork, info);
	 }
	 template< typename T, typename Ta >
	 void pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, T* ab, integer* ldab, T* afb, 
               integer* ldafb, char* equed, Ta* s, T* b, integer* ldb, T* x, integer* ldx, Ta* rcond, Ta* ferr,
               Ta* berr, T* work, Ta* rwork, integer* info)
	 {
	   pbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work,
            rwork, info);
	 }

	 template< typename T >
	 void pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  T* ab, integer* ldab, T* b, integer* ldb,
               integer* info)
	 {
	   pbtrs(uplo, n, kd, nrhs,  ab, ldab, b, ldb, info);
	 }

	 template< typename T >
	 void pocon(char* uplo, integer* n,  T* a, integer* lda, T* anorm, T* rcond, T* work, integer* iwork,
               integer* info)
	 {
	   pocon(uplo, n,  a, lda, anorm, rcond, work, iwork, info);
	 }
	 template< typename T, typename Ta >
	 void pocon(char* uplo, integer* n,  T* a, integer* lda, Ta* anorm, Ta* rcond, T* work, Ta* rwork, 
               integer* info)
	 {
	   pocon(uplo, n,  a, lda, anorm, rcond, work, rwork, info);
	 }

	 template< typename T >
	 void posv(char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* b, integer* ldb, integer* info)
	 {
	   posv(uplo, n, nrhs, a, lda, b, ldb, info);
	 }

	 template< typename T >
	 void posvx(char* fact, char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* af, integer* ldaf,
               char* equed, T* s, T* b, integer* ldb, T* x, integer* ldx, T* rcond, T* ferr, T* berr, T* work,
               integer* iwork, integer* info)
	 {
	   posvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork,
            info);
	 }
	 template< typename T, typename Ta >
	 void posvx(char* fact, char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* af, integer* ldaf,
               char* equed, Ta* s, T* b, integer* ldb, T* x, integer* ldx, Ta* rcond, Ta* ferr, Ta* berr,
               T* work, Ta* rwork, integer* info)
	 {
	   posvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork,
            info);
	 }

	 template< typename T >
	 void posvxx(char* fact, char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* af, integer* ldaf,
                char* equed, T* s, T* b, integer* ldb, T* x, integer* ldx, T* rcond, T* rpvgrw, T* berr, 
                integer* n_err_bnds, T* err_bnds_norm, T* err_bnds_comp, integer* nparams, T* params, 
                float* work, integer* iwork, integer* info)
	 {
	   posvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds,
             err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
	 }
	 template< typename T, typename Ta >
	 void posvxx(char* fact, char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* af, integer* ldaf,
                char* equed, Ta* s, T* b, integer* ldb, T* x, integer* ldx, Ta* rcond, Ta* rpvgrw, Ta* berr,
                integer* n_err_bnds, Ta* err_bnds_norm, Ta* err_bnds_comp, integer* nparams, Ta* params,
                dcomplex* work, double* rwork, integer* info)
	 {
	   posvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds,
             err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
	 }

	 template< typename T >
	 void ppsv(char* uplo, integer* n, integer* nrhs, T* ap, T* b, integer* ldb, integer* info)
	 {
	   ppsv(uplo, n, nrhs, ap, b, ldb, info);
	 }

	 template< typename T >
	 void ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, T* ap, T* afp, char* equed, T* s, T* b,
               integer* ldb, T* x, integer* ldx, T* rcond, T* ferr, T* berr, T* work, integer* iwork, 
               integer* info)
	 {
	   ppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
	 }
	 template< typename T, typename Ta >
	 void ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, T* ap, T* afp, char* equed, Ta* s, T* b,
               integer* ldb, T* x, integer* ldx, Ta* rcond, Ta* ferr, Ta* berr, T* work, Ta* rwork, 
               integer* info)
	 {
	   ppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
	 }

	 template< typename T >
	 void ptsv(integer* n, integer* nrhs, T* d, T* e, T* b, integer* ldb, integer* info)
	 {
	   ptsv(n, nrhs, d, e, b, ldb, info);
	 }
	 template< typename T, typename Ta >
	 void ptsv(integer* n, integer* nrhs, Ta* d, T* e, T* b, integer* ldb, integer* info)
	 {
	   ptsv(n, nrhs, d, e, b, ldb, info);
	 }

	 template< typename T >
	 void ptsvx(char* fact, integer* n, integer* nrhs, T* d,  T* e, T* df, T* ef, T* b, integer* ldb,
               T* x, integer* ldx, T* rcond, T* ferr, T* berr, T* work, integer* info)
	 {
	   ptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, info);
	 }
	 template< typename T, typename Ta >
	 void ptsvx(char* fact, integer* n, integer* nrhs, Ta* d, T* e, Ta* df, T* ef, T* b, integer* ldb,
               T* x, integer* ldx, Ta* rcond, Ta* ferr, Ta* berr, T* work, Ta* rwork, integer* info)
	 {
	   ptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
	 }

	 template< typename T >
	 void sbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* w, T* z,
                     integer* ldz, T* work, integer* lwork, integer* info)
	 {
	   sbev_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info);
	 }
	 template< typename T, typename Ta >
	 void hbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, Ta* w, T* z,
                     integer* ldz, T* work, integer* lwork, Ta* rwork, integer* info)
	 {
	   hbev_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, info);
	 }

	 template< typename T >
	 void sbev(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* w, T* z, integer* ldz,
              T* work, integer* info)
	 {
	   sbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info);
	 }
	 template< typename T, typename Ta >
	 void hbev(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, Ta* w, T* z, integer* ldz, 
              T* work, Ta* rwork, integer* info)
	 {
	   hbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info);
	 }

	 template< typename T >
	 void sbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* w, T* z, 
                      integer* ldz, T* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
	 {
	   sbevd_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
	 }
	 template< typename T, typename Ta >
	 void hbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, Ta* w, T* z, 
                      integer* ldz, T* work, integer* lwork, Ta* rwork, integer* lrwork, integer* iwork, 
                      integer* liwork, integer* info)
	 {
	   hbevd_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
	 }

	 template< typename T >
	 void sbevd(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* w, T* z, integer* ldz, 
               T* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
	 {
	   sbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
	 }
	 template< typename T, typename Ta >
	 void hbevd(char* jobz, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, Ta* w, T* z, integer* ldz, 
               T* work, integer* lwork, Ta* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
	 {
	   hbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
	 }

	 template< typename T >
	 void sbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* q,
                      integer* ldq, T* vl, T* vu, integer* il, integer* iu, T* abstol, integer* m, T* w, T* z,
                      integer* ldz, T* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
	 {
	   sbevx_2stage(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                   iwork, ifail, info);
	 }
	 template< typename T, typename Ta >
	 void hbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* q,
                      integer* ldq, Ta* vl, Ta* vu, integer* il, integer* iu, Ta* abstol, integer* m, Ta* w, T* z,
                      integer* ldz, T* work, integer* lwork, Ta* rwork, integer* iwork, integer* ifail, integer* info)
	 {
	   hbevx_2stage(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                   rwork, iwork, ifail, info);
	 }

	 template< typename T >
	 void sbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* q,
               integer* ldq, T* vl, T* vu, integer* il, integer* iu, T* abstol, integer* m, T* w, T* z,
               integer* ldz, T* work, integer* iwork, integer* ifail, integer* info)
	 {
	   sbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, 
            info);
	 }
	 template< typename T, typename Ta >
	 void hbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, T* ab, integer* ldab, T* q, 
               integer* ldq, Ta* vl, Ta* vu, integer* il, integer* iu, Ta* abstol, integer* m, Ta* w, T* z, 
               integer* ldz, T* work, Ta* rwork, integer* iwork, integer* ifail, integer* info)
	 {
	   hbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork,
            ifail, info);
	 }

	 template< typename T >
	 void sbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, T* ab, integer* ldab,  T* bb, 
               integer* ldbb, T* x, integer* ldx, T* work, integer* info)
	 {
	   sbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, info);
	 }
	 template< typename T, typename Ta >
	 void hbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, T* ab, integer* ldab,  T* bb,
               integer* ldbb, T* x, integer* ldx, T* work, Ta* rwork, integer* info)
	 {
	   hbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info);
	 }

	 template< typename T >
	 void syevd_2stage(char* jobz, char* uplo, integer* n, T* a, integer* lda, T* w, T* work, integer* lwork,
                      integer* iwork, integer* liwork, integer* info)
	 {
	   syevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
	 }

	 template< typename T >
	 void sytrs_rook(char* uplo, integer* n, integer* nrhs,  T* a, integer* lda, integer* ipiv, T* b, 
                    integer* ldb, integer* info)
	 {
	   sytrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
	 }

	 template< typename T, typename Ta >
	 void hecon_3(char* uplo, integer* n, T* a, integer* lda, T* e, integer* ipiv, Ta* anorm, Ta* rcond,
                 T* work, integer* info)
	 {
	   hecon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info);
	 }

	 template< typename T, typename Ta >
	 void heequb(char* uplo, integer* n, T* a, integer* lda, Ta* s, Ta* scond, Ta* amax, T* work, integer* info)
	 {
	   heequb(uplo, n, a, lda, s, scond, amax, work, info);
	 }

	 template< typename T, typename Ta >
	 void heev_2stage(char* jobz, char* uplo, integer* n, T* a, integer* lda, Ta* w, T* work, integer* lwork,
                     Ta* rwork, integer* info)
	 {
	   heev_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
	 }

	 template< typename T, typename Ta >
	 void heevd_2stage(char* jobz, char* uplo, integer* n, T* a, integer* lda, Ta* w, T* work, integer* lwork, 
                      Ta* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
	 {
	   heevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
	 }

	 template< typename T, typename Ta >
	 void heevr_2stage(char* jobz, char* range, char* uplo, integer* n, T* a, integer* lda, Ta* vl, Ta* vu,
                      integer* il, integer* iu, Ta* abstol, integer* m, Ta* w, T* z, integer* ldz, 
                      integer* isuppz, T* work, integer* lwork, Ta* rwork, integer* lrwork, integer* iwork, 
                      integer* liwork, integer* info)
	 {
	   heevr_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, 
                   rwork, lrwork, iwork, liwork, info);
	 }

	 template< typename T, typename Ta >
	 void heevx(char* jobz, char* range, char* uplo, integer* n, T* a, integer* lda, Ta* vl, Ta* vu, integer* il,
               integer* iu, Ta* abstol, integer* m, Ta* w, T* z, integer* ldz, T* work, integer* lwork, 
               Ta* rwork, integer* iwork, integer* ifail, integer* info)
	 {
	   heevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, 
            ifail, info);
	 }

	 template< typename T, typename Ta >
	 void heevx_2stage(char* jobz, char* range, char* uplo, integer* n, T* a, integer* lda, Ta* vl, Ta* vu, 
                      integer* il, integer* iu, Ta* abstol, integer* m, Ta* w, T* z, integer* ldz, T* work, 
                      integer* lwork, Ta* rwork, integer* iwork, integer* ifail, integer* info)
	 {
	   heevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, 
                   iwork, ifail, info);
	 }

	 template< typename T, typename Ta >
	 void hegv(integer* itype, char* jobz, char* uplo, integer* n, T* a, integer* lda, T* b, integer* ldb, 
              Ta* w, T* work, integer* lwork, Ta* rwork, integer* info)
	 {
	   hegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
	 }

	 template< typename T, typename Ta >
	 void hegv_2stage(integer* itype, char* jobz, char* uplo, integer* n, T* a, integer* lda, T* b,
                     integer* ldb, Ta* w, T* work, integer* lwork, Ta* rwork, integer* info)
	 {
	   hegv_2stage(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
	 }

	 template< typename T, typename Ta >
	 void hegvd(integer* itype, char* jobz, char* uplo, integer* n, T* a, integer* lda, T* b, integer* ldb,
               Ta* w, T* work, integer* lwork, Ta* rwork, integer* lrwork, integer* iwork, integer* liwork, 
               integer* info)
	 {
	   hegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info);
	 }

	 template< typename T, typename Ta >
	 void hegvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, T* a, integer* lda, T* b, 
               integer* ldb, Ta* vl, Ta* vu, integer* il, integer* iu, Ta* abstol, integer* m, Ta* w, T* z,
               integer* ldz, T* work, integer* lwork, Ta* rwork, integer* iwork, integer* ifail, integer* info)
	 {
	   hegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, 
            rwork, iwork, ifail, info);
	 }

	 template< typename T, typename Ta >
	 void herfs(char* uplo, integer* n, integer* nrhs,  T* a, integer* lda,  T* af, integer* ldaf, integer* ipiv,
               T* b, integer* ldb, T* x, integer* ldx, Ta* ferr, Ta* berr, T* work, Ta* rwork, integer* info)
	 {
	   herfs(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
	 }

	 template< typename T, typename Ta >
	 void herfsx(char* uplo, char* equed, integer* n, integer* nrhs, T* a, integer* lda, T* af, integer* ldaf, 
                integer* ipiv, Ta* s, T* b, integer* ldb, T* x, integer* ldx, Ta* rcond, Ta* berr, 
                integer* n_err_bnds, Ta* err_bnds_norm, Ta* err_bnds_comp, integer* nparams, Ta* params, 
                T* work, Ta* rwork, integer* info)
	 {
	   herfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, 
             err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
	 }

	 template< typename T >
	 void hesv(char* uplo, integer* n, integer* nrhs, T* a, integer* lda, integer* ipiv, T* b, integer* ldb,
              T* work, integer* lwork, integer* info)
	 {
	   hesv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
	 }

	 template< typename T >
	 void hesv_aa(char* uplo, integer* n, integer* nrhs, T* a, integer* lda, integer* ipiv, T* b, integer* ldb,
                 T* work, integer* lwork, integer* info)
	 {
	   hesv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
	 }

	 template< typename T >
	 void hesv_aa_2stage(char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* tb, integer* ltb,
                        integer* ipiv, integer* ipiv2, T* b, integer* ldb, T* work, integer* lwork, 
                        integer* info)
	 {
	   hesv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
	 }

	 template< typename T >
	 void hesv_rk(char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* e, integer* ipiv, T* b, 
                 integer* ldb, T* work, integer* lwork, integer* info)
	 {
	   hesv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
	 }

	 template< typename T, typename Ta >
	 void hesvx(char* fact, char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* af, integer* ldaf, 
               integer* ipiv, T* b, integer* ldb, T* x, integer* ldx, Ta* rcond, Ta* ferr, Ta* berr, T* work, 
               integer* lwork, Ta* rwork, integer* info)
	 {
	   hesvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork,
            info);
	 }

	 template< typename T, typename Ta >
	 void hesvxx(char* fact, char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* af, integer* ldaf, 
                integer* ipiv, char* equed, Ta* s, T* b, integer* ldb, T* x, integer* ldx, Ta* rcond, 
                Ta* rpvgrw, Ta* berr, integer* n_err_bnds, Ta* err_bnds_norm, Ta* err_bnds_comp, 
                integer* nparams, Ta* params, T* work, Ta* rwork, integer* info)
	 {
	   hesvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, 
             n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
	 }

	 template< typename T >
	 void heswapr(char* uplo, integer* n, T* a, integer* lda, integer* i1, integer* i2)
	 {
	   heswapr(uplo, n, a, lda, i1, i2);
	 }

    template< typename T >
	 void hetrf(char* uplo, integer* n, T* a, integer* lda, integer* ipiv, T* work, integer* lwork, 
               integer* info)
	 {
	   hetrf(uplo, n, a, lda, ipiv, work, lwork, info);
	 }

	 template< typename T >
	 void hetrf_rk(char* uplo, integer* n, T* a, integer* lda, T* e, integer* ipiv, T* work, integer* lwork, 
                  integer* info)
	 {
	   hetrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info);
	 }

	 template< typename T >
	 void hetri(char* uplo, integer* n, T* a, integer* lda, integer* ipiv, T* work, integer* info)
	 {
	   hetri(uplo, n, a, lda, ipiv, work, info);
	 }

	 template< typename T >
	 void hetri_3(char* uplo, integer* n, T* a, integer* lda,  T* e,  integer* ipiv, T* work, integer* lwork, 
                 integer* info)
	 {
	   hetri_3(uplo, n, a, lda, e, ipiv, work, lwork, info);
	 }

	 template< typename T >
	 void hetri2(char* uplo, integer* n, T* a, integer* lda, integer* ipiv, T* work, integer* lwork,
                integer* info)
	 {
	   hetri2(uplo, n, a, lda, ipiv, work, lwork, info);
	 }

	 template< typename T >
	 void hetri2x(char* uplo, integer* n, T* a, integer* lda, integer* ipiv, T* work, integer *nb, integer* info)
	 {
	   hetri2x(uplo, n, a, lda, ipiv, work, nb, info);
	 }

	 template< typename T >
	 void hetrs(char* uplo, integer* n, integer* nrhs,  T* a, integer* lda,  integer* ipiv, T* b, integer* ldb,
               integer* info)
	 {
	   hetrs(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
	 }

	 template< typename T >
	 void hetrs_3(char* uplo, integer* n, integer* nrhs, T* a, integer* lda, T* e, integer* ipiv, T* b, 
                 integer* ldb, integer* info)
	 {
	   hetrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info);
	 }

	 template< typename T >
	 void hetrs_rook(char* uplo, integer* n, integer* nrhs,  T* a, integer* lda,  integer* ipiv, T* b, 
                    integer* ldb, integer* info)
	 {
	   hetrs_rook(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
	 }

    template< typename T >
	 void hetrs2(char* uplo, integer* n, integer* nrhs,  T* a, integer* lda,  integer* ipiv, T* b, integer* ldb,
                T* work, integer* info)
	 {
	   hetrs2(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, info);
	 }

	 template< typename T >
	 void combssq(T* v1, T* v2)
	 {
	   combssq(v1, v2);
	 }

	 template< typename T >
	 void hesv_rook(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b, 
                   integer *ldb, T *work, integer *lwork, integer *info)
	 {
	   hesv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
	 }

	 template< typename T >
	 void hetf2(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, integer *info)
	 {
	   hetf2(uplo, n, a, lda, ipiv, info); 
	 }

	 template< typename T >
	 void hetf2_rk(char *uplo, integer *n, T *a, integer *lda, T *e, integer *ipiv, integer *info)
	 {
	   hetf2_rk(uplo, n, a, lda, e, ipiv, info); 
	 }

	 template< typename T >
	 void hetf2_rook(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, integer *info)
	 {
	   hetf2_rook(uplo, n, a, lda, ipiv, info); 
	 }

	 template< typename T, typename Ta >
	 void hetrd_hb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, T *ab, integer *ldab,
                     Ta *d, Ta * e, T *hous, integer *lhous, T *work, integer *lwork, integer *info)
	 {
	   hetrd_hb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
	 }

	 template< typename T >
	 void hetri_3x(char *uplo, integer *n, T *a, integer * lda, T *e, integer *ipiv, T *work, integer *nb, 
                  integer * info)
	 {
	   hetri_3x(uplo, n, a, lda, e, ipiv, work, nb, info);
	 }

	 template< typename T >
	 void hetri_rook(char *uplo, integer *n, T *a, integer * lda, integer *ipiv, T *work, integer * info)
	 {
	   hetri_rook(uplo, n, a, lda, ipiv, work, info);
	 }

	 template< typename T >
	 void pbtf2(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, integer *info)
	 {
	   pbtf2(uplo, n, kd, ab, ldab, info);
	 }

	 template< typename T >
	 void spmv(char *uplo, integer *n, T *alpha, T *ap, T *x, integer *incx, T *beta, T *y, integer * incy)
	 {
	   spmv(uplo, n, alpha, ap, x, incx, beta, y, incy);
	 }

	 template< typename T >
	 void spr(char *uplo, integer *n, T *alpha, T *x, integer *incx, T *ap)
	 {
	   spr(uplo, n, alpha, x, incx, ap);
	 }

	 template< typename T >
	 void symv(char *uplo, integer *n, T *alpha, T * a, integer *lda, T *x, integer *incx, T *beta, T *y,
              integer *incy)
	 {
	   symv(uplo, n, alpha, a, lda, x, incx, beta, y, incy); 
	 }

	 template< typename T >
	 void sytrd_sy2sb(char *uplo, integer *n, integer *kd, T *a, integer *lda, T *ab, integer *ldab, T *tau,
                     T *work, integer *lwork, integer *info)
	 {
	   sytrd_sy2sb(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
	 }

	 template< typename T >
	 void la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, 
                            T *af, integer * ldaf, logical *colequ, T *c, T *b, integer *ldb, T *y, 
                            integer *ldy, T *berr_out, integer *n_norms, T * err_bnds_norm, T *err_bnds_comp, 
                            T *res, T *ayb, T * dy, T *y_tail, T *rcond, integer *ithresh, T *rthresh, T *dz_ub, 
                            logical *ignore_cwise, integer *info)
	 {
	   la_porfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                         err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, 
                         ignore_cwise, info);
	 }
	 template< typename T, typename Ta >
	 void la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                            integer *ldaf, logical *colequ, Ta *c, T *b, integer *ldb, T *y, integer *ldy, 
                            Ta *berr_out, integer *n_norms, Ta * err_bnds_norm, Ta *err_bnds_comp, T *res, 
                            Ta *ayb, T *dy, T *y_tail, Ta *rcond, integer *ithresh, Ta * rthresh, Ta *dz_ub, 
                            logical *ignore_cwise, integer *info)
	 {
	   la_porfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, 
                         n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, 
                         dz_ub, ignore_cwise, info);
	 }

	 template< typename T >
	 void la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, 
                            T *af, integer * ldaf, integer *ipiv, logical *colequ, T *c, T *b, integer * ldb, 
                            T *y, integer *ldy, T *berr_out, integer *n_norms, T *err_bnds_norm, 
                            T *err_bnds_comp, T *res, T *ayb, T *dy, T *y_tail, T *rcond, integer *ithresh, 
                            T * rthresh, T *dz_ub, logical *ignore_cwise, integer *info)
	 {
	   la_syrfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out,
                         n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, 
                         dz_ub, ignore_cwise, info);
	 }
	 template< typename T, typename Ta >
	 void la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, 
                            T *af, integer *ldaf, integer *ipiv, logical *colequ, Ta *c, T *b, integer *ldb, 
                            T *y, integer *ldy, Ta *berr_out, integer * n_norms, Ta *err_bnds_norm, 
                            Ta *err_bnds_comp, T *res, Ta *ayb, T *dy, T *y_tail, Ta *rcond, integer * ithresh, 
                            Ta *rthresh, Ta *dz_ub, logical *ignore_cwise, integer *info)
	 {
	   la_syrfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, 
                         n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, 
                         dz_ub, ignore_cwise, info);
	 }

	 template< typename T >
	 void la_syamv(integer *uplo, integer *n, T *alpha, T *a, integer *lda, T *x, integer *incx, T *beta, T *y, 
                  integer *incy)
	 {
	   la_syamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
	 }
	 template< typename T, typename Ta >
	 void la_syamv(integer *uplo, integer *n, Ta *alpha, T *a, integer *lda, T *x, integer *incx, Ta *beta, 
                  Ta *y, integer *incy)
	 {
	   la_syamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
	 }

	 inline void cladiv(scomplex *ret_val, scomplex *x, scomplex *y)
	 {
	   cladiv_(ret_val, x, y);
	 }
	 inline void zladiv(dcomplex *ret_val, dcomplex *x, dcomplex *y)
	 {
	   zladiv_(ret_val, x, y);
	 }

	 template< typename T >
	 void lahef(char *uplo, integer *n, integer *nb, integer *kb, T *a, integer *lda, integer *ipiv, T *w,
               integer *ldw, integer *info)
	 {
	   lahef(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
	 }

	 template< typename T >
	 void lahef_rk(char *uplo, integer *n, integer *nb, integer *kb, T *a, integer *lda, T *e, integer *ipiv,
                  T *w, integer *ldw, integer *info)
	 {
	   lahef_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
	 }

	 template< typename T >
	 void lahef_rook(char *uplo, integer *n, integer *nb, integer *kb, T *a, integer *lda, integer *ipiv, T *w, 
                    integer *ldw, integer *info)
	 {
	   lahef_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info); 
	 }

	 template< typename T >
	 void lahrd(integer *n, integer *k, integer *nb, T *a, integer *lda, T *tau, T *t, integer *ldt, T *y, 
               integer *ldy)
	 {
	   lahrd(n, k, nb, a, lda, tau, t, ldt, y, ldy);
	 }

	 template< typename T >
	 void lamrg(integer *n1, integer *n2, T *a, integer * strd1, integer *strd2, integer *index)
	 {
	   lamrg(n1, n2, a, strd1, strd2, index);
	 }

	 template< typename T >
	 void laqsb(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, T *s, T *scond, T *amax, char *equed)
	 {
	   laqsb(uplo, n, kd, ab, ldab, s, scond, amax, equed);
	 }
	 template< typename T, typename Ta >
	 void laqsb(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, Ta *s, Ta *scond, Ta *amax, 
               char *equed)
	 {
	   laqsb(uplo, n, kd, ab, ldab, s, scond, amax, equed);
	 }

	 template< typename T >
	 void laqsp(char *uplo, integer *n, T *ap, T *s, T * scond, T *amax, char *equed)
	 {
	   laqsp(uplo, n, ap, s, scond, amax, equed);
	 }
	 template< typename T, typename Ta >
	 void laqsp(char *uplo, integer *n, T *ap, Ta *s, Ta * scond, Ta *amax, char *equed)
	 {
	   laqsp(uplo, n, ap, s, scond, amax, equed);
	 }

	 template< typename T >
	 void laqsy(char *uplo, integer *n, T *a, integer *lda, T *s, T *scond, T *amax, char *equed)
	 {
	   laqsy(uplo, n, a, lda, s, scond, amax, equed);
	 }
	 template< typename T, typename Ta >
	 void laqsy(char *uplo, integer *n, T *a, integer *lda, Ta *s, Ta *scond, Ta *amax, char *equed)
	 {
	   laqsy(uplo, n, a, lda, s, scond, amax, equed);
	 }

	 template< typename T >
	 void larrr(integer *n, T *d, T *e, integer *info)
	 {
	   larrr(n, d, e, info);   
	 }

	 template< typename T >
	 void lasd5(integer *i, T *d, T *z, T *delta, T *rho, T *dsigma, T *work)
	 {
	   lasd5(i, d, z, delta, rho, dsigma, work); 
	 }

	 template< typename T >
	 void laswlq(integer *m, integer *n, integer *mb, integer * nb, T *a, integer *lda, T *t, integer *ldt, 
                T *work, integer *lwork, integer *info)
	 {
	   laswlq(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
	 }

	 template< typename T >
	 void lasy2(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2, T *tl, integer *ldtl, 
               T *tr, integer * ldtr, T *b, integer *ldb, T *scale, T *x, integer *ldx, T *xnorm, integer *info)
	 {
	   lasy2(ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr,  ldtr, b, ldb, scale, x, ldx, xnorm, info);
	 }

	 inline void dlat2s(char *uplo, integer *n, double *a, integer * lda, float *sa, integer *ldsa, integer *info)
	 {
	   dlat2s_(uplo, n, a, lda, sa, ldsa, info); 
	 }

	 inline void zlat2c(char *uplo, integer *n, dcomplex *a, integer *lda, scomplex *sa, integer *ldsa, 
                       integer *info)
	 {
	   zlat2c_(uplo, n, a, lda, sa, ldsa, info);
	 }

	 template< typename T >
	 void gegs(char* jobvsl, char* jobvsr, integer* n, T* a, integer* lda, T* b, integer* ldb, T* alphar, 
              T* alphai, T* beta, T* vsl, integer* ldvsl, T* vsr, integer* ldvsr, T* work, integer* lwork,
              integer* info)
	 {
	   gegs(jobvsl, jobvsr, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, info);
	 }
	 template <typename T, typename Ta>
	 void gegs(char* jobvsl, char* jobvsr, integer* n, T* a, integer* lda, T* b, integer* ldb, T* alpha, T* beta,
              T* vsl, integer* ldvsl, T* vsr, integer* ldvsr, T* work, integer* lwork, Ta* rwork, integer* info)
	 {
	   gegs(jobvsl, jobvsr, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, info);
	 }

	 template< typename T >
	 void gegv(char* jobvl, char* jobvr, integer* n, T* a, integer* lda, T* b, integer* ldb, T* alphar, 
              T* alphai, T* beta, T* vl, integer* ldvl, T* vr, integer* ldvr, T* work, integer* lwork, 
              integer* info)
	 {
	   ggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
	 }
	 template <typename T, typename Ta>
	 void gegv(char* jobvl, char* jobvr, integer* n, T* a, integer* lda, T* b, integer* ldb, T* alpha, T* beta, 
              T* vl, integer* ldvl, T* vr, integer* ldvr, T* work, integer* lwork, Ta* rwork, integer* info)
	 {
	   gegv(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
	 }

	 template< typename T >
	 void latzm(char* side, integer* m, integer* n, T* v, integer* incv, T* tau, T* c1, T* c2, integer* ldc, 
               T* work)
	 {
	   latzm(side, m, n, v, incv, tau, c1, c2, ldc, work);
	 }

	 inline void dsposv(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, 
                       double* x, integer* ldx, double* work, float* swork, integer* iter, integer* info)
	 {
	   dsposv_(uplo, n, nrhs, a, lda, b, ldb, x, ldx, work, swork, iter, info);
	 }

	 inline void zcposv(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, 
                       integer* ldb, dcomplex* x, integer* ldx, dcomplex* work, scomplex* swork, double* rwork, 
                       integer* iter, integer* info)
	 {
	   zcposv_(uplo, n, nrhs, a, lda, b, ldb, x, ldx, work, swork, rwork, iter, info);
	 }

	 template< typename T >
	 void getrfnp(integer* m, integer* n, T* a, integer* lda, integer* info)
	 {
	   getrfnp(m, n, a, lda, info);
	 }

	 template< typename T >
	 void getrfnpi(integer *m, integer *n, integer *nfact, T *a, integer *lda, integer *info)
	 {
	   getrfnpi(m, n, nfact, a, lda, info);
	 }

	 template< typename T >
	 void spffrt2(T *ap, integer *n, integer * ncolm, T *work, T *work2)
	 {
	   spffrt2(ap, n, ncolm, work, work2);
	 }
	 template< typename T >
	 void spffrtx(T *ap, integer *n, integer * ncolm, T *work, T *work2)
	 {
	   spffrtx(ap, n, ncolm, work, work2);
	 }
}
#endif