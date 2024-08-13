/******************************************************************************
* Copyright (C) 2021-2024, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame.hh
 *  libflame.hh defines all the overloaded CPP functions to be invoked from
 *  template interfaces
 *  */
#ifndef LIBFLAME_HH
#define LIBFLAME_HH

#include  <FLAME.h>

namespace libflame {

// Cholesky factorization of a real symmetric positive definite matrix a
inline void potrf(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  spotrf_(uplo, n, a, lda, info);
}
inline void potrf(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  dpotrf_(uplo, n, a, lda, info);
}
inline void potrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  cpotrf_(uplo, n, a, lda, info);
}
inline void potrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  zpotrf_(uplo, n, a, lda, info);
}

// --- computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (unblocked algorithm) ---
inline void potf2(char* uplo, integer*n, float* a, integer* lda, integer* info)
{
  spotf2_(uplo, n, a, lda, info);
}
inline void potf2(char* uplo, integer*n, double* a, integer* lda, integer* info)
{
  dpotf2_(uplo, n, a, lda, info);
}
inline void potf2(char* uplo, integer*n, scomplex* a, integer* lda, integer* info)
{
  cpotf2_(uplo, n, a, lda, info);
}
inline void potf2(char* uplo, integer*n, dcomplex* a, integer* lda, integer* info)
{
  zpotf2_(uplo, n, a, lda, info);
}

// --- LU factorization with partial pivoting
inline void getrf(integer* m, integer* n, float* a, integer* lda, integer* ipiv, integer* info)
{
  sgetrf_(m, n, a, lda, ipiv, info);
}
inline void getrf(integer* m, integer* n, double*   a, integer* lda, integer* ipiv, integer* info)
{
  dgetrf_(m, n, a, lda, ipiv, info);
}
inline void getrf(integer* m, integer* n, scomplex* a, integer* lda, integer* ipiv, integer* info)
{
  cgetrf_(m, n, a, lda, ipiv, info);
}
inline void getrf(integer* m, integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* info)
{
  zgetrf_(m, n, a, lda, ipiv, info);
}


//--- computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm) ---
inline void getf2(integer* m, integer* n, float* a, integer* lda, integer* ipiv, integer* info)
{
  sgetf2_(m, n, a, lda, ipiv, info);
}
inline void getf2(integer* m, integer* n, double* a, integer* lda, integer* ipiv, integer* info)
{
  dgetf2_(m, n, a, lda, ipiv, info);
}
inline void getf2(integer* m, integer* n, scomplex* a, integer* lda, integer* ipiv, integer* info)
{
  cgetf2_(m, n, a, lda, ipiv, info);
}
inline void getf2(integer* m, integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* info)
{
  zgetf2_(m, n, a, lda, ipiv, info);
}

// --- QR factorization (classic) ---
inline void geqrf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqrf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqrf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqrf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
  
//--- computes the QR factorization of a general rectangular matrix using an unblocked algorithm ---  
inline void geqr2(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* info)
{
  sgeqr2_(m, n, a, lda, tau, work, info);
}
inline void geqr2(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* info)
{
  dgeqr2_(m, n, a, lda, tau, work, info);
}
inline void geqr2(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* info)
{
  cgeqr2_(m, n, a, lda, tau, work, info);
}
inline void geqr2(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* info)
{
  zgeqr2_(m, n, a, lda, tau, work, info);
}
  
//--- QR factorization with column pivoting of a real m-by-n matrix a ---  
inline void geqpf(integer* m, integer* n, float* a, integer* lda, integer* jpvt, float* tau, float* work, integer* info)
{
  printf(" Function sgeqpf() has been deprecated. Please use sgep3() instead.\n");
  sgeqpf_(m, n, a, lda, jpvt, tau, work, info);
}
inline void geqpf(integer* m, integer* n, double* a, integer* lda, integer* jpvt, double* tau, double* work, integer* info)
{
  printf(" Function dgeqpf() has been deprecated. Please use dgep3() instead.\n");
  dgeqpf_(m, n, a, lda, jpvt, tau, work, info);
}
inline void geqpf(integer* m, integer* n, scomplex* a, integer* lda, integer* jpvt, scomplex* tau, scomplex* work, float* rwork, integer* info)
{
  printf(" Function cgeqpf() has been deprecated. Please use cgep3() instead.\n");
  cgeqpf_(m, n, a, lda, jpvt, tau, work, rwork, info);
}
inline void geqpf(integer* m, integer* n, dcomplex* a, integer* lda, integer* jpvt, dcomplex* tau, dcomplex* work, double* rwork, integer* info)
{
  printf(" Function zgeqpf() has been deprecated. Please use zgep3() instead.\n");
  zgeqpf_(m, n, a, lda, jpvt, tau, work, rwork, info);
}

// --- QR factorization with column pivoting of a real m-by-n matrix ---
inline void geqp3(integer* m, integer* n, float* a, integer* lda, integer* jpvt, float* tau, float* work, integer* lwork, integer* info)
{
  sgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
}
inline void geqp3(integer* m, integer* n, double* a, integer* lda, integer* jpvt, double* tau, double* work, integer* lwork, integer* info)
{
  dgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
}
inline void geqp3(integer* m, integer* n, scomplex* a, integer* lda, integer* jpvt, scomplex* tau, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
}
inline void geqp3(integer* m, integer* n, dcomplex* a, integer* lda, integer* jpvt, dcomplex* tau, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
}

// --- LQ factorization (classic) ---
inline void gelqf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sgelqf_(m, n, a, lda, tau, work, lwork, info);
}
inline void gelqf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dgelqf_(m, n, a, lda, tau, work, lwork, info);
}
inline void gelqf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cgelqf_(m, n, a, lda, tau, work, lwork, info);
}
inline void gelqf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zgelqf_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes the LQ factorization of a general rectangular matrix using an unblocked algorithm ---
inline void gelq2(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* info)
{
  sgelq2_(m, n, a, lda, tau, work, info);
}
inline void gelq2(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* info)
{
  dgelq2_(m, n, a, lda, tau, work, info);
}
inline void gelq2(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* info)
{
  cgelq2_(m, n, a, lda, tau, work, info);
}
inline void gelq2(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* info)
{
  zgelq2_(m, n, a, lda, tau, work, info);
}

// --- LS solver ---
inline void gelsd(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float*  s, float*  rcond, integer* rank, float* work, integer* lwork, integer* iwork, integer* info)
{
  sgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
}
inline void gelsd(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* s, double* rcond, integer* rank, double* work, integer* lwork, integer* iwork, integer* info)
{
  dgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
}
inline void gelsd(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, float*  s, float*  rcond, integer* rank, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  cgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);
}
inline void gelsd(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* s, double* rcond, integer* rank, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  zgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);
}

// ---  minimum-norm solution to a real linear least squares problem ---
inline void gelss(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float*  s, float*  rcond, integer* rank, float* work, integer* lwork, integer* info)
{
  sgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}
inline void gelss(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* s, double* rcond, integer* rank, double* work, integer* lwork, integer* info)
{
  dgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}
inline void gelss(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, float*  s, float*  rcond, integer* rank, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);
}
inline void gelss(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* s, double* rcond, integer* rank, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);
}

// --- Triangular-transpose matrix multiply ---
inline void lauum(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  slauum_(uplo, n, a, lda, info);
}
inline void lauum(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  dlauum_(uplo, n, a, lda, info);
}
inline void lauum(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  clauum_(uplo, n, a, lda, info);
}
inline void lauum(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  zlauum_(uplo, n, a, lda, info);
}

//--- computes the product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm) ---
inline void lauu2(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  slauu2_(uplo, n, a, lda, info);
}
inline void lauu2(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  dlauu2_(uplo, n, a, lda, info);
}
inline void lauu2(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  clauu2_(uplo, n, a, lda, info);
}
inline void lauu2(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  zlauu2_(uplo, n, a, lda, info);
}
// --- Symmetric (hermitian) positive definite matrix inversion ---
inline void potri(char* uplo, integer*  n, float* buff_A, integer*  ldim_A, integer* info)
{
  spotri_(uplo, n, buff_A, ldim_A, info);
}
inline void potri(char* uplo, integer*  n, double* buff_A, integer*  ldim_A, integer* info)
{
  dpotri_(uplo, n, buff_A, ldim_A, info);
}
inline void potri(char* uplo, integer*  n, scomplex* buff_A, integer*  ldim_A, integer* info)
{
  cpotri_(uplo, n, buff_A, ldim_A, info);
}
inline void potri(char* uplo, integer*  n, dcomplex* buff_A, integer*  ldim_A, integer* info)
{
  zpotri_(uplo, n, buff_A, ldim_A, info);
}

// --- Triangular matrix inversion ---
inline void trtri(char* uplo, char* diag, integer* n, float* a, integer* lda, integer* info)
{
  strtri_(uplo, diag, n, a, lda, info);
}
inline void trtri(char* uplo, char* diag, integer* n, double* a, integer* lda, integer* info)
{
  dtrtri_(uplo, diag, n, a, lda, info);
}
inline void trtri(char* uplo, char* diag, integer* n, scomplex* a, integer* lda, integer* info)
{
  ctrtri_(uplo, diag, n, a, lda, info);
}
inline void trtri(char* uplo, char* diag, integer* n, dcomplex* a, integer* lda, integer* info)
{
  ztrtri_(uplo, diag, n, a, lda, info);
}

//--- computes the inverse of a triangular matrix (unblocked algorithm)---
inline void trti2(char* uplo, char* diag, integer* n, float* a, integer* lda, integer* info)
{
  strti2_(uplo, diag, n, a, lda, info);
}
inline void trti2(char* uplo, char* diag, integer* n, double* a, integer* lda, integer* info)
{
  dtrti2_(uplo, diag, n, a, lda, info);
}
inline void trti2(char* uplo, char* diag, integer* n, scomplex* a, integer* lda, integer* info)
{
  ctrti2_(uplo, diag, n, a, lda, info);
}
inline void trti2(char* uplo, char* diag, integer* n, dcomplex* a, integer* lda, integer* info)
{
  ztrti2_(uplo, diag, n, a, lda, info);
}

// --- Triangular Sylvester equation solve ---
inline void trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, float* a, integer* lda, float* b, integer* ldb, float* c, integer* ldc, float* scale, integer* info)
{
  strsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}
inline void trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, double* a, integer* lda, double* b, integer* ldb, double* c, integer* ldc, double* scale, integer* info)
{
  dtrsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}
inline void trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* c, integer* ldc, float* scale, integer* info)
{
  ctrsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}
inline void trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* c, integer* ldc, double* scale, integer* info)
{
  ztrsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}

// --- Reduction to upper Hessenberg form ---
inline void gehrd(integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}
inline void gehrd(integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}
inline void gehrd(integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}
inline void gehrd(integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}

// --- reduces a general square matrix to upper Hessenberg form using an unblocked algorithm ---
inline void gehd2(integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* tau, float* work, integer* info)
{
  sgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}
inline void gehd2(integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* tau, double* work, integer* info)
{
  dgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}
inline void gehd2(integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* info)
{
  cgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}
inline void gehd2(integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* info)
{
  zgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}

// --- Reduction to tridiagonal form ---
inline void sytrd(char* uplo, integer* n, float* a, integer* lda, float*  d, float*  e, float* tau, float* work, integer* lwork, integer* info)
{
  ssytrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}
inline void sytrd(char* uplo, integer* n, double* a, integer* lda, double* d, double* e, double* tau, double* work, integer* lwork, integer* info)
{
  dsytrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}

//---reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation ---
inline void sytd2(char* uplo, integer* n, float* a, integer* lda, float*  d, float*  e, float* tau, integer* info)
{
  ssytd2_(uplo, n, a, lda, d, e, tau, info);
}
inline void sytd2(char* uplo, integer* n, double* a, integer* lda, double* d, double* e, double* tau, integer* info)
{
  dsytd2_(uplo, n, a, lda, d, e, tau, info);
}

//--- reduces a Hermitian matrix to real symmetric tridiagonal form by an unitary similarity transformation ---
inline void hetd2(char* uplo, integer* n, scomplex* a, integer* lda, float*  d, float*  e, scomplex* tau, integer* info)
{
  chetd2_(uplo, n, a, lda, d, e, tau, info);
}
inline void hetd2(char* uplo, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tau, integer* info)
{
  zhetd2_(uplo, n, a, lda, d, e, tau, info);
}

// --- Reduction to bidiagonal form ---
inline void gebrd(integer* m, integer* n, float* a, integer* lda, float*  d, float*  e, float* tauq, float* taup, float* work, integer* lwork, integer* info)
{
  sgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}
inline void gebrd(integer* m, integer* n, double* a, integer* lda, double* d, double* e, double* tauq, double* taup, double* work, integer* lwork, integer* info)
{
  dgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}
inline void gebrd(integer* m, integer* n, scomplex* a, integer* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, integer* lwork, integer* info)
{
  cgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}
inline void gebrd(integer* m, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, integer* lwork, integer* info)
{
  zgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}

// --- reduces a general matrix to bidiagonal form using an unblocked algorithm ---
inline void gebd2(integer* m, integer* n, float* a, integer* lda, float*  d, float*  e, float* tauq, float* taup, float* work, integer* info)
{
  sgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}
inline void gebd2(integer* m, integer* n, double* a, integer* lda, double* d, double* e, double* tauq, double* taup, double* work, integer* info)
{
  dgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}
inline void gebd2(integer* m, integer* n, scomplex* a, integer* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, integer* info)
{
  cgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}
inline void gebd2(integer* m, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, integer* info)
{
  zgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}

// --- Reduce Hermitian-definite generalized eigenproblem to standard form ---
inline void sygst(integer* itype, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  ssygst_(itype, uplo, n, a, lda, b, ldb, info);
}
inline void sygst(integer* itype, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  dsygst_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- reduces a complex Hermitian-definite generalized eigenproblem to standard form ---
inline void hegst(integer* itype, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  chegst_(itype, uplo, n, a, lda, b, ldb, info);
}
inline void hegst(integer* itype, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  zhegst_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- reduces a symmetric definite generalized eigenproblem to standard form ---
inline void sygs2(integer* itype, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  ssygs2_(itype, uplo, n, a, lda, b, ldb, info);
}
inline void sygs2(integer* itype, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  dsygs2_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- reduces a Hermitian definite generalized eigenproblem to standard form ---
inline void hegs2(integer* itype, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  chegs2_(itype, uplo, n, a, lda, b, ldb, info);
}
inline void hegs2(integer* itype, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  zhegs2_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- Accumulate block Householder matrix T (classic) ---
inline void larft(char* direct, char* storev, integer* n, integer* k, float* v, integer* ldv, float* tau, float* t, integer* ldt)
{
  slarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}
inline void larft(char* direct, char* storev, integer* n, integer* k, double* v, integer* ldv, double* tau, double* t, integer* ldt)
{
  dlarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}
inline void larft(char* direct, char* storev, integer* n, integer* k, scomplex* v, integer* ldv, scomplex* tau, scomplex* t, integer* ldt)
{
  clarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}
inline void larft(char* direct, char* storev, integer* n, integer* k, dcomplex* v, integer* ldv, dcomplex* tau, dcomplex* t, integer* ldt)
{
  zlarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}

// --- Generate a Householder vector (classic) ---
inline void larfg( integer* n, float* alpha, float* x, integer* incx, float* tau)
{
  slarfg_(n, alpha, x, incx, tau);
}
inline void larfg( integer* n, double* alpha, double* x, integer* incx, double* tau)
{
  dlarfg_(n, alpha, x, incx, tau);
}
inline void larfg( integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* tau)
{
  clarfg_(n, alpha, x, incx, tau);
}
inline void larfg( integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* tau)
{
  zlarfg_(n, alpha, x, incx, tau);
}

// --- generates an elementary reflector (Householder matrix) with non-negative beta ---
inline void larfgp( integer* n, float* alpha, float* x, integer* incx, float* tau)
{
  slarfgp_(n, alpha, x, incx, tau);
}
inline void larfgp( integer* n, double* alpha, double* x, integer* incx, double* tau)
{
  dlarfgp_(n, alpha, x, incx, tau);
}
inline void larfgp( integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* tau)
{
  clarfgp_(n, alpha, x, incx, tau);
}
inline void larfgp( integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* tau)
{
  zlarfgp_(n, alpha, x, incx, tau);
}

// --- Form Q from QR factorization ---
inline void orgqr(integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sorgqr_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void orgqr(integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dorgqr_(m, n, k, a, lda, tau, work, lwork, info);
}
// --- generates an M-by-N complex matrix Q with orthonormal columns ---
inline void ungqr(integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cungqr_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void ungqr(integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zungqr_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from QR factorization ---
inline void ormqr(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline void ormqr(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

inline void unmqr(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmqr_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmqr(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// ---multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf ---  
inline void orm2r(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* info)
{
  sorm2r_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void orm2r(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* info)
{
  dorm2r_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

inline void unm2r(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  cunm2r_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, info);
}
inline void unm2r(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  zunm2r_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- Form Q from LQ factorization ---
inline void orglq(integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sorglq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void orglq(integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dorglq_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- generates an M-by-N complex matrix Q ---
inline void unglq(integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cunglq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void unglq(integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zunglq_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from LQ factorization ---
inline void ormlq(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sormlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline void ormlq(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dormlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- overwrites the general complex M-by-N matrix C ---
inline void unmlq(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmlq(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a LQ factorization determined by sgelqf ---
inline void orml2(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* info)
{
  sorml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void orml2(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* info)
{
  dorml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

inline void unml2(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  cunml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void unml2(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  zunml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- Form Q from tridiagonal reduction ---
inline void orgtr(char* uplo, integer* m, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sorgtr_(uplo, m, a, lda, tau, work, lwork, info);
}
inline void orgtr(char* uplo, integer* m, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dorgtr_(uplo, m, a, lda, tau, work, lwork, info);
}

// --- generates a complex unitary matrix Q ---
inline void ungtr(char* uplo, integer* m, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cungtr_(uplo, m, a, lda, tau, work, lwork, info);
}
inline void ungtr(char* uplo, integer* m, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zungtr_(uplo, m, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from tridiagonal reduction ---
inline void ormtr(char* side, char* uplo, char* trans, integer* m, integer* n, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sormtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}
inline void ormtr(char* side, char* uplo, char* trans, integer* m, integer* n, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dormtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}

// --- UNMTR overwrites the general complex M-by-N matrix C ---
inline void unmtr(char* side, char* uplo, char* trans, integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmtr(char* side, char* uplo, char* trans, integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}

// --- Form Q from bidiagonal reduction ---
inline void orgbr(char* vect, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sorgbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}
inline void orgbr(char* vect, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dorgbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}

// --- generates one of the complex unitary matrices Q ---  
inline void ungbr(char* vect, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cungbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}
inline void ungbr(char* vect, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zungbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from bidiagonal reduction ---
inline void ormbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sormbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline void ormbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dormbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- overwrites the general complex M-by-N matrix C ---
inline void unmbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- Tridiagonal QR algorithm ---
inline void steqr(char* compz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* info)
{
  ssteqr_(compz, n, d, e, z, ldz, work, info);
}
inline void steqr(char* compz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* info)
{
  dsteqr_(compz, n, d, e, z, ldz, work, info);
}
inline void steqr(char* compz, integer* n, float* d, float* e, scomplex* z, integer* ldz, float* work, integer* info)
{
  csteqr_(compz, n, d, e, z, ldz, work, info);
}
inline void steqr(char* compz, integer* n, double* d, double* e, dcomplex* z, integer* ldz, double* work, integer* info)
{
  zsteqr_(compz, n, d, e, z, ldz, work, info);
}

// --- Tridiagonal divide-and-conquer algorithm ---
inline void stedc(char* compz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  sstedc_(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}
inline void stedc(char* compz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dstedc_(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}
inline void stedc(char* compz, integer* n, float* d, float* e, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  cstedc_(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void stedc(char* compz, integer* n, double* d, double* e, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zstedc_(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- Tridiagonal MRRR algorithm ---
inline void stemr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, integer* m, float* w, float* z, integer* ldz, integer* nzc, integer* isuppz, logical* tryrac, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  sstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}
inline void stemr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, integer* m, double* w, double* z, integer* ldz, integer* nzc, integer* isuppz, integer* tryrac, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}
inline void stemr(char* jobz, char* range, integer* n, float*  d, float*  e, float* vl, float* vu, integer* il, integer* iu, integer* m, float*  w, scomplex* z, integer* ldz, integer* nzc, integer* isuppz, integer* tryrac, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  cstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}
inline void stemr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, integer* m, double* w, dcomplex* z, integer* ldz, integer* nzc, integer* isuppz, integer* tryrac, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  zstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}

// --- Hermitian eigenvalue decomposition (QR algorithm) ---
inline void syev(char* jobz, char* uplo, integer* n, float* a, integer* lda, float*  w, float* work, integer* lwork, integer* info)
{
  ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}
inline void syev(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, integer* info)
{
  dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

// --- performs the symmetric rank-1 update of a complex symmetric matrix ---
inline void syr(char* uplo, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* a, integer* lda)
{
  csyr_(uplo, n, alpha, x, incx, a, lda);
}
inline void syr(char* uplo, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* a, integer* lda)
{
  zsyr_(uplo, n, alpha, x, incx, a, lda);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heev(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}
inline void heev(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heevd(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float*  w, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  cheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void heevd(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heevr(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float*  vl, float*  vu, integer* il, integer* iu, float*  abstol, integer* m, float*  w, scomplex* z, integer* ldz, integer* isuppz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  cheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void heevr(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, integer* isuppz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- Hermitian eigenvalue decomposition (divide-and-conquer) ---
inline void syevd(char* jobz, char* uplo, integer* n, float* a, integer* lda, float*  w, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}
inline void syevd(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}

// --- Hermitian eigenvalue decomposition (MRRR) ---
inline void syevr(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssyevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline void syevr(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsyevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- Bidiagonal QR algorithm ---
inline void bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, float* d, float* e, float* vt, integer* ldvt, float* u, integer* ldu, float* c, integer* ldc, float* rwork, integer* info)
{
  sbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}
inline void bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, double* d, double* e, double* vt, integer* ldvt, double* u, integer* ldu, double* c, integer* ldc, double* rwork, integer* info)
{
  dbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}
inline void bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, float* d, float* e, scomplex* vt, integer* ldvt, scomplex* u, integer* ldu, scomplex* c, integer* ldc, float* rwork, integer* info)
{
  cbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}
inline void bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, double* d, double* e, dcomplex* vt, integer* ldvt, dcomplex* u, integer* ldu, dcomplex* c, integer* ldc, double* rwork, integer* info)
{
  zbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}

// --- Bidiagonal divide-and-conquor algorithm ---
inline void bdsdc(char* uplo, char* compq, integer* n, float*  d, float*  e, float*  u, integer* ldu, float*  vt, integer* ldvt, float*  q, float*  iq, float* work, integer* iwork, integer* info)
{
  sbdsdc_(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info);
}
inline void bdsdc(char* uplo, char* compq, integer* n, double* d, double* e, double* u, integer* ldu, double* vt, integer* ldvt, double* q, double* iq, double* work, integer* iwork, integer* info)
{
  dbdsdc_(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info);
}

// --- computes the singular value decomposition (SVD) of a real N-by-N (upper or lower) bidiagonal matrix ---
inline void bdsvdx(char* uplo, char* jobz, char* range, integer* n, float* d, float* e, float *vl, float *vu, integer* il, integer* iu, integer* ns, float* s, float* z, integer* ldz, float* work, integer* iwork, integer* info)
{
  sbdsvdx_(uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work, iwork, info);
}
inline void bdsvdx(char* uplo, char* jobz, char* range, integer* n, double* d, double* e, double *vl, double *vu, integer* il, integer* iu, integer* ns, double* s, double* z, integer* ldz, double* work, integer* iwork, integer* info)
{
  dbdsvdx_(uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work, iwork, info);
}

// --- computes the reciprocal condition numbers for the eigenvectors of a real symmetric or complex Hermitian matrix ---
inline void disna(char* job, integer* m, integer* n,  float* d, float* sep, integer* info)
{
  sdisna_(job, m, n,  d, sep, info);
}
inline void disna(char* job, integer* m, integer* n,  double* d, double* sep, integer* info)
{
  ddisna_(job, m, n,  d, sep, info);
}

// --- General matrix singular value decomposition (QR algorithm) ---
inline void gesvd(char* jobu, char* jobv, integer* m, integer* n, float* a, integer* lda, float* s, float* u, integer* ldu, float* vt, integer* ldvt, float* work, integer* lwork, integer* info)
{
  sgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}
inline void gesvd(char* jobu, char* jobv, integer* m, integer* n, double* a, integer* lda, double* s, double* u, integer* ldu, double* vt, integer* ldvt, double* work, integer* lwork, integer* info)
{
  dgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}
inline void gesvd(char* jobu, char* jobv, integer* m, integer* n, scomplex* a, integer* lda, float*  s, scomplex* u, integer* ldu, scomplex* vt, integer* ldvt, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}
inline void gesvd(char* jobu, char* jobv, integer* m, integer* n, dcomplex* a, integer* lda, double* s, dcomplex* u, integer* ldu, dcomplex* vt, integer* ldvt, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}

// --- General matrix singular value decomposition (divide-and-conquer) ---
inline void gesdd(char* jobz, integer* m, integer* n, float* a, integer* lda, float*  s, float* u, integer* ldu, float* vt, integer* ldvt, float* work, integer* lwork, integer* iwork, integer* info)
{
  sgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline void gesdd(char* jobz, integer* m, integer* n, double* a, integer* lda, double* s, double* u, integer* ldu, double* vt, integer* ldvt, double* work, integer* lwork, integer* iwork, integer* info)
{
  dgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline void gesdd(char* jobz, integer* m, integer* n, scomplex* a, integer* lda, float*  s, scomplex* u, integer* ldu, scomplex* vt, integer* ldvt, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  cgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}
inline void gesdd(char* jobz, integer* m, integer* n, dcomplex* a, integer* lda, double* s, dcomplex* u, integer* ldu, dcomplex* vt, integer* ldvt, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  zgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}

// --- Swap rows ---
inline void laswp(integer* n, float* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  slaswp_(n, a, lda, k1, k2, ipiv, incx);
}
inline void laswp(integer* n, double* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  dlaswp_(n, a, lda, k1, k2, ipiv, incx);
}
inline void laswp(integer* n, scomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  claswp_(n, a, lda, k1, k2, ipiv, incx);
}
inline void laswp(integer* n, dcomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  zlaswp_(n, a, lda, k1, k2, ipiv, incx);
}

// --- Initialize a matrix ---
inline void laset(char* uplo, integer* m, integer* n, float* alpha, float* beta, float* a, integer* lda)
{
  slaset_(uplo, m, n, alpha, beta, a, lda);
}
inline void laset(char* uplo, integer* m, integer* n, double* alpha, double* beta, double* a, integer* lda)
{
  dlaset_(uplo, m, n, alpha, beta, a, lda);
}
inline void laset(char* uplo, integer* m, integer* n, scomplex* alpha, scomplex* beta, scomplex* a, integer* lda)
{
  claset_(uplo, m, n, alpha, beta, a, lda);
}
inline void laset(char* uplo, integer* m, integer* n, dcomplex* alpha, dcomplex* beta, dcomplex* a, integer* lda)
{
  zlaset_(uplo, m, n, alpha, beta, a, lda);
}

  // --- Bidiagonal block cs decomposition of orthogonal/unitary matrix  ---
inline void bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, float* theta, float* phi, float* u1, integer* ldu1, float* u2, integer* ldu2, float* v1t, integer* ldv1t, float* v2t, integer* ldv2t, float* b11d, float* b11e, float* b12d, float* b12e, float* b21d, float* b21e, float* b22d, float* b22e, float* rwork, integer* lrwork, integer* info)
{
  sbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
inline void bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, double* theta, double* phi, double* u1, integer* ldu1, double* u2, integer* ldu2, double* v1t, integer* ldv1t, double* v2t, integer* ldv2t, double* b11d, double* b11e, double* b12d, double* b12e, double* b21d, double* b21e, double* b22d, double* b22e, double* rwork, integer* lrwork, integer* info)
{
  dbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
inline void bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, float* theta, float* phi, scomplex* u1, integer* ldu1, scomplex* u2, integer* ldu2, scomplex* v1t, integer* ldv1t, scomplex* v2t, integer* ldv2t, float* b11d, float* b11e, float* b12d, float* b12e, float* b21d, float* b21e, float* b22d, float* b22e, float* rwork, integer* lrwork, integer* info)
{
  cbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
inline void bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, double* theta, double* phi, dcomplex* u1, integer* ldu1, dcomplex* u2, integer* ldu2, dcomplex* v1t, integer* ldv1t, dcomplex* v2t, integer* ldv2t, double* b11d, double* b11e, double* b12d, double* b12e, double* b21d, double* b21e, double* b22d, double* b22e, double* rwork, integer* lrwork, integer* info)
{
  zbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
  // ---  Reduces a general band matrix to bidiagonal form.  ---
inline void gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, float* ab, integer* ldab, float* d, float* e, float* q, integer* ldq, float* pt, integer* ldpt, float* c, integer* ldc, float* work, integer* info)
{
  sgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
}
inline void gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, double* ab, integer* ldab, double* d, double* e, double* q, integer* ldq, double* pt, integer* ldpt, double* c, integer* ldc, double* work, integer* info)
{
  dgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
}
inline void gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, scomplex* ab, integer* ldab, float* d, float* e, scomplex* q, integer* ldq, scomplex* pt, integer* ldpt, scomplex* c, integer* ldc, scomplex* work, float* rwork, integer* info)
{
  cgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, rwork, info);
}
inline void gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, dcomplex* ab, integer* ldab, double* d, double* e, dcomplex* q, integer* ldq, dcomplex* pt, integer* ldpt, dcomplex* c, integer* ldc, dcomplex* work, double* rwork, integer* info)
{
  zgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, rwork, info);
}

// ---estimates the reciprocal of the condition number of a real general band matrix A, in either the 1-norm or the infinity-norm ---
inline void gbcon(char* norm, integer* n, integer* kl, integer* ku, float* ab, integer* ldab, integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  sgbcon_(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info);
}
inline void gbcon(char* norm, integer* n, integer* kl, integer* ku, double* ab, integer* ldab,  integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dgbcon_(norm, n, kl, ku, ab, ldab,  ipiv, anorm, rcond, work, iwork, info);
}
inline void gbcon(char* norm, integer* n, integer* kl, integer* ku, scomplex* ab, integer* ldab,  integer* ipiv, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  cgbcon_(norm, n, kl, ku, ab, ldab,  ipiv, anorm, rcond, work, rwork, info);
}
inline void gbcon(char* norm, integer* n, integer* kl, integer* ku, dcomplex* ab, integer* ldab,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  zgbcon_(norm, n, kl, ku, ab, ldab,  ipiv, anorm, rcond, work, rwork, info);
}

// --- Computes row and column scalings intended to equilibrate  band matrix and reduce its condition number---
inline void gbequ(integer* m, integer* n, integer* kl, integer* ku,  float* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  sgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline void gbequ(integer* m, integer* n, integer* kl, integer* ku,  double* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  dgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline void gbequ(integer* m, integer* n, integer* kl, integer* ku,  scomplex* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  cgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline void gbequ(integer* m, integer* n, integer* kl, integer* ku,  dcomplex* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  zgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}

// --- Computes row and column scalings intended to equilibrate  band matrix and reduce its condition number(scaling factor:power of radix)---
inline void gbequb(integer* m, integer* n, integer* kl, integer* ku,  float* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  sgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline void gbequb(integer* m, integer* n, integer* kl, integer* ku,  double* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  dgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline void gbequb(integer* m, integer* n, integer* kl, integer* ku,  scomplex* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  cgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline void gbequb(integer* m, integer* n, integer* kl, integer* ku,  dcomplex* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  zgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}

// ---Improves computed solution to banded matrix and provides error bounds and backward estimates ---
inline void gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, integer* ipiv, float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  double* ab, integer* ldab, double* afb, integer* ldafb,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// ---Improves computed solution to banded matrix and provides backward estimates,component wise and normwise error bounds---
inline void gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  float* ab, integer* ldab,  float* afb, integer* ldafb,  integer* ipiv,  float* r,  float* c,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  sgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  double* ab, integer* ldab,  double* afb, integer* ldafb,  integer* ipiv,  double* r,  double* c,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  scomplex* ab, integer* ldab,  scomplex* afb, integer* ldafb,  integer* ipiv,  float* r,  float* c,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  cgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  dcomplex* ab, integer* ldab,  dcomplex* afb, integer* ldafb,  integer* ipiv,  double* r,  double* c,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// ---Computes the solution to system of linear equations A * X = B for GB matrices ---
inline void gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, integer* ipiv, float* b, integer* ldb, integer* info)
{
  sgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
inline void gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, double* ab, integer* ldab, integer* ipiv, double* b, integer* ldb, integer* info)
{
  dgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
inline void gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, scomplex* ab, integer* ldab, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  cgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
inline void gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, dcomplex* ab, integer* ldab, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}

// ---It computes the solution to system of linear equations A * X = B for GB matrices along with error bounds ---
inline void gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, double* ab, integer* ldab, double* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// ---computes the solution to system of linear equations A * X = B for GB matrices with error, normalisations ---
inline void gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  sgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, 
  rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, double* ab, integer* ldab, double* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
  ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  cgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, 
  rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// ---Computes the LU factorization of a general band matrix using partial pivot---
inline void gbtrf(integer* m, integer* n, integer* kl, integer* ku, float* ab, integer* ldab, integer* ipiv, integer* info)
{
  sgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}
inline void gbtrf(integer* m, integer* n, integer* kl, integer* ku, double* ab, integer* ldab, integer* ipiv, integer* info)
{
  dgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}
inline void gbtrf(integer* m, integer* n, integer* kl, integer* ku, scomplex* ab, integer* ldab, integer* ipiv, integer* info)
{
  cgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}
inline void gbtrf(integer* m, integer* n, integer* kl, integer* ku, dcomplex* ab, integer* ldab, integer* ipiv, integer* info)
{
  zgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}

// ---solves a system of linear equationsA * X = B  or  A**T * X = B with a general band matrix A using the LU factorization computed by GBTRF ---
inline void gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  float* ab, integer* ldab,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  sgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}
inline void gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  double* ab, integer* ldab,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  dgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}
inline void gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  scomplex* ab, integer* ldab,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  cgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}
inline void gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  dcomplex* ab, integer* ldab,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}

// ---Forms right left eigen vectors of real general matrix ---
inline void gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* scale, integer* m, float* v, integer* ldv, integer* info)
{
  sgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}
inline void gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* scale, integer* m, double* v, integer* ldv, integer* info)
{
  dgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}
inline void gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* scale, integer* m, scomplex* v, integer* ldv, integer* info)
{
  cgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}
inline void gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* scale, integer* m, dcomplex* v, integer* ldv, integer* info)
{
  zgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}

// --- Balances a real matrix ---
inline void gebal(char* job, integer* n, float* a, integer* lda, integer* ilo, integer* ihi, float* scale, integer* info)
{
  sgebal_(job, n, a, lda, ilo, ihi, scale, info);
}
inline void gebal(char* job, integer* n, double* a, integer* lda, integer* ilo, integer* ihi, double* scale, integer* info)
{
  dgebal_(job, n, a, lda, ilo, ihi, scale, info);
}
inline void gebal(char* job, integer* n, scomplex* a, integer* lda, integer* ilo, integer* ihi, float* scale, integer* info)
{
  cgebal_(job, n, a, lda, ilo, ihi, scale, info);
}
inline void gebal(char* job, integer* n, dcomplex* a, integer* lda, integer* ilo, integer* ihi, double* scale, integer* info)
{
  zgebal_(job, n, a, lda, ilo, ihi, scale, info);
}

// ---Computes row and column scaling to reduce condition number of matrix ---
inline void geequb(integer* m, integer* n,  float* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  sgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline void geequb(integer* m, integer* n,  double* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  dgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline void geequb(integer* m, integer* n,  scomplex* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  cgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline void geequb(integer* m, integer* n,  dcomplex* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  zgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}

// ---Computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices ---
inline void gees(char* jobvs, char* sort, void* select, integer* n, float* a, integer* lda, integer* sdim, float* wr, float* wi, float* vs, integer* ldvs, float* work, integer* lwork, logical* bwork, integer* info)
{
  sgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info);
}
inline void gees(char* jobvs, char* sort, void* select, integer* n, double* a, integer* lda, integer* sdim, double* wr, double* wi, double* vs, integer* ldvs, double* work, integer* lwork, logical* bwork, integer* info)
{
  dgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info);
}
inline void gees(char* jobvs, char* sort, void* select, integer* n, scomplex* a, integer* lda, integer* sdim, scomplex* w, scomplex* vs, integer* ldvs, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  cgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info);
}
inline void gees(char* jobvs, char* sort, void* select, integer* n, dcomplex* a, integer* lda, integer* sdim, dcomplex* w, dcomplex* vs, integer* ldvs, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  zgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info);
}

// ---computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices with rconde and rcondv---
inline void geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, float* a, integer* lda, integer* sdim, float* wr, float* wi, float* vs, integer* ldvs, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  sgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline void geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, double* a, integer* lda, integer* sdim, double* wr, double* wi, double* vs, integer* ldvs, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  dgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline void geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, scomplex* a, integer* lda, integer* sdim, scomplex* w, scomplex* vs, integer* ldvs, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  cgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, w, vs, ldvs, rconde, rcondv, work, lwork, rwork, bwork, info);
}
inline void geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, dcomplex* a, integer* lda, integer* sdim, dcomplex* w, dcomplex* vs, integer* ldvs, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  zgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, w, vs, ldvs, rconde, rcondv, work, lwork, rwork, bwork, info);
}

// ---Computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices ---
inline void geev(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* wr, float* wi, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  sgeev_(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void geev(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* wr, double* wi, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  dgeev_(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void geev(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* w, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cgeev_(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline void geev(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* w, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zgeev_(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices(enabling conditions) ---
inline void geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, float* a, integer* lda, float* wr, float* wi, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* ilo, integer* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, integer* info)
{
  sgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
}
inline void geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, double* a, integer* lda, double* wr, double* wi, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* ilo, integer* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, integer* info)
{
  dgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
}
inline void geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, scomplex* a, integer* lda, scomplex* w, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* ilo, integer* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
}
inline void geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, dcomplex* a, integer* lda, dcomplex* w, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* ilo, integer* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
}

// ---Computes the singular value decomposition (SVD) of a real matrix ---
inline void gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
             char* jobp, integer* m, integer* n, float* a, integer* lda, float* sva,
             float* u, integer* ldu, float* v, integer* ldv,
             float* work, integer* lwork, integer* iwork, integer* info)
{
  sgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork, info);
}
inline void gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
             char* jobp, integer* m, integer* n, double* a, integer* lda, double* sva,
             double* u, integer* ldu, double* v, integer* ldv, double* work, integer* lwork, integer* iwork, integer* info)
{
  dgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv,  work, lwork, iwork, info);
}
inline void gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, integer* m,
             integer* n, scomplex* a, integer* lda, float* sva, scomplex* u, integer* ldu,
             scomplex* v, integer* ldv, scomplex* cwork, integer* lwork, float* rwork, integer* lrwork, 
             integer* iwork, integer* info)
{
  cgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, cwork, lwork, rwork, lrwork, iwork, info);
}
inline void gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, integer* m, integer* n, dcomplex* a, integer* lda, double* sva, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* cwork, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* info)
{
  zgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, cwork, lwork, rwork, lrwork, iwork, info);
}

// --- Computes LQ factorization of a real matrix ---
inline void gelq(integer* m, integer* n, float* a, integer* lda, float* t, integer* tsize, float* work, integer* lwork, integer* info)
{
  sgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline void gelq(integer* m, integer* n, double* a, integer* lda, double* t, integer* tsize, double* work, integer* lwork, integer* info)
{
  dgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline void gelq(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* tsize, scomplex* work, integer* lwork, integer* info)
{
  cgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline void gelq(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* tsize, dcomplex* work, integer* lwork, integer* info)
{
  zgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}

// --- Solves overdetermined or underdetermined systems for GE matrices ---
inline void gels(char* trans, integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  sgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline void gels(char* trans, integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline void gels(char* trans, integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  cgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline void gels(char* trans, integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}

// --- Solves overdetermined or underdetermined systems for GE matrices ---
inline void gelsy(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, float* work, integer* lwork, integer* info)
{
  sgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
}
inline void gelsy(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, double* work, integer* lwork, integer* info)
{
  dgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
}
inline void gelsy(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);
}
inline void gelsy(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);
}

inline void gelsx(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, float* work, integer* info)
{
  printf(" Function sgelsx() has been deprecated. Please use sgelsy() instead.\n"); 
  sgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, info);
}
inline void gelsx(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, double* work, integer* info)
{
  printf(" Function dgelsx() has been deprecated. Please use dgelsy() instead.\n"); 
  dgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, info);
}
inline void gelsx(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, scomplex* work, float* rwork, integer* info)
{
  printf(" Function cgelsx() has been deprecated. Please use cgelsy() instead.\n"); 
  cgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, rwork, info);
}
inline void gelsx(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, dcomplex* work, double* rwork, integer* info)
{
  printf(" Function zgelsx() has been deprecated. Please use zgelsy() instead.\n"); 
  zgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, rwork, info);
}

// --- Overwrites general matrix with a form compatible with orthogonal matrix ---
inline void gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  float* a, integer* lda,  float* t, integer* tsize, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline void gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  double* a, integer* lda,  double* t, integer* tsize, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline void gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  scomplex* a, integer* lda,  scomplex* t, integer* tsize, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline void gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  dcomplex* a, integer* lda,  dcomplex* t, integer* tsize, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}

// --- Multiples a matrix C by a real orthogonal or complex unitary matrix Q, as computed by ?geqr---
inline void gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  float* a, integer* lda,  float* t, integer* tsize, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline void gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  double* a, integer* lda,  double* t, integer* tsize, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline void gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  scomplex* a, integer* lda,  scomplex* t, integer* tsize, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline void gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  dcomplex* a, integer* lda,  dcomplex* t, integer* tsize, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}

// ---Multiplies a general matrix by the orthogonal/unitary matrix Q of the QR factorization formed by geqrt. ---
inline void gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  float* v, integer* ldv,  float* t, integer* ldt, float* c, integer* ldc, float* work, integer* info)
{
  sgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}
inline void gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  double* v, integer* ldv,  double* t, integer* ldt, double* c, integer* ldc, double* work, integer* info)
{
  dgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}
inline void gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  cgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}
inline void gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  zgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}

// --- computes a QL factorization of M-by-N matrix ---
inline void geqlf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sgeqlf_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqlf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dgeqlf_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqlf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cgeqlf_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqlf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zgeqlf_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes a QR factorization of M-by-N matrix ---
inline void geqr(integer* m, integer* n, float* a, integer* lda, float* t, integer* tsize, float* work, integer* lwork, integer* info)
{
  sgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline void geqr(integer* m, integer* n, double* a, integer* lda, double* t, integer* tsize, double* work, integer* lwork, integer* info)
{
  dgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline void geqr(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* tsize, scomplex* work, integer* lwork, integer* info)
{
  cgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline void geqr(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* tsize, dcomplex* work, integer* lwork, integer* info)
{
  zgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}

// --- computes a QR factorization of a M-by-N matrix ---
inline void geqrfp(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqrfp(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqrfp(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}
inline void geqrfp(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes a blocked QR factorization---
inline void geqrt(integer* m, integer* n, integer* nb, float* a, integer* lda, float* t, integer* ldt, float* work, integer* info)
{
  sgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}
inline void geqrt(integer* m, integer* n, integer* nb, double* a, integer* lda, double* t, integer* ldt, double* work, integer* info)
{
  dgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}
inline void geqrt(integer* m, integer* n, integer* nb, scomplex* a, integer* lda, scomplex* t, integer* ldt, scomplex* work, integer* info)
{
  cgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}
inline void geqrt(integer* m, integer* n, integer* nb, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, dcomplex* work, integer* info)
{
  zgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}

// --- computes a blocked QR factorization ---
inline void geqrt2(integer* m, integer* n, float* a, integer* lda, float* t, integer* ldt, integer* info)
{
  sgeqrt2_(m, n, a, lda, t, ldt, info);
}
inline void geqrt2(integer* m, integer* n, double* a, integer* lda, double* t, integer* ldt, integer* info)
{
  dgeqrt2_(m, n, a, lda, t, ldt, info);
}
inline void geqrt2(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* ldt, integer* info)
{
  cgeqrt2_(m, n, a, lda, t, ldt, info);
}
inline void geqrt2(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, integer* info)
{
  zgeqrt2_(m, n, a, lda, t, ldt, info);
}

// --- computes a blocked QR factorization ---
inline void geqrt3(integer* m, integer* n, float* a, integer* lda, float* t, integer* ldt, integer* info)
{
  sgeqrt3_(m, n, a, lda, t, ldt, info);
}
inline void geqrt3(integer* m, integer* n, double* a, integer* lda, double* t, integer* ldt, integer* info)
{
  dgeqrt3_(m, n, a, lda, t, ldt, info);
}
inline void geqrt3(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* ldt, integer* info)
{
  cgeqrt3_(m, n, a, lda, t, ldt, info);
}
inline void geqrt3(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, integer* info)
{
  zgeqrt3_(m, n, a, lda, t, ldt, info);
}

// --- Improves the computed solution to a system of linear equations---
inline void gerfs(char* trans, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void gerfs(char* trans, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void gerfs(char* trans, integer* n, integer* nrhs, scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void gerfs(char* trans, integer* n, integer* nrhs, dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- GERFSX improves the computed solution to a system of linear equations ---
inline void gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* r,  float* c,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  sgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* r,  double* c,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  float* r,  float* c,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  cgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  double* r,  double* c,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes a RQ factorization of a M-by-N matrix ---
inline void gerqf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sgerqf_(m, n, a, lda, tau, work, lwork, info);
}
inline void gerqf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dgerqf_(m, n, a, lda, tau, work, lwork, info);
}
inline void gerqf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cgerqf_(m, n, a, lda, tau, work, lwork, info);
}
inline void gerqf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zgerqf_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes the solution to a real system of linear equations ---
inline void gesv(integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, integer* info)
{
  sgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline void gesv(integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, integer* info)
{
  dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline void gesv(integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  cgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline void gesv(integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

// --- computes the singular value decomposition (SVD) ---
inline void gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, integer* n, float* a, integer* lda, float* s, float* u, integer* ldu, float* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork, float* work, integer* lwork, float* rwork, integer* lrwork, integer* info)
{
  sgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, 
            iwork, liwork, work, lwork, rwork, lrwork, info);
}
inline void gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, integer* n, double* a, integer* lda, double* s, double* u, integer* ldu, double* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork, double* work, integer* lwork, double* rwork, integer* lrwork, integer* info)
{
  dgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork, liwork, work, lwork, rwork, lrwork, info);
}
inline void gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, 
            integer* n, scomplex* a, integer* lda, float* s, scomplex* u, integer* ldu,
            scomplex* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork,
            scomplex* cwork, integer* lcwork, float* rwork, integer* lrwork, integer* info)
{
  cgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork, liwork, cwork, lcwork, rwork, lrwork, info);
}
inline void gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, integer* n, dcomplex* a, integer* lda, double* s, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork, dcomplex* cwork, integer* lcwork, double* rwork, integer* lrwork, integer* info)
{
  zgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork, liwork, cwork, lcwork, rwork, lrwork, info);
}

// --- computes the singular value decomposition (SVD) for GE matrices ---
inline void gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, integer* ns, float* s, float* u, integer* ldu, float* vt, integer* ldvt, float* work, integer* lwork, integer* iwork, integer* info)
{
  sgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline void gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, integer* ns, double* s, double* u, integer* ldu, double* vt, integer* ldvt, double* work, integer* lwork, integer* iwork, integer* info)
{
  dgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline void gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, integer* ns, float* s, scomplex* u, integer* ldu, scomplex* vt, integer* ldvt, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  cgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}
inline void gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, integer* ns, double* s, dcomplex* u, integer* ldu, dcomplex* vt, integer* ldvt, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  zgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}

// --- computes the singular value decomposition (SVD) of a M-by-N matrix ---
inline void gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, float* a, integer* lda, float* sva, integer* mv, float* v, integer* ldv, float* work, integer* lwork, integer* info)
{
  sgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info);
}
inline void gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, double* a, integer* lda, double* sva, integer* mv, double* v, integer* ldv, double* work, integer* lwork, integer* info)
{
  dgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info);
}
inline void gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, scomplex* a, integer* lda, float* sva, integer* mv, scomplex* v, integer* ldv, scomplex* cwork, integer* lwork, float* rwork, integer* lrwork, integer* info)
{
  cgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork, lwork, rwork, lrwork, info);
}
inline void gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, dcomplex* a, integer* lda, double* sva, integer* mv, dcomplex* v, integer* ldv, dcomplex* cwork, integer* lwork, double* rwork, integer* lrwork, integer* info)
{
  zgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork, lwork, rwork, lrwork, info);
}

// --- computes the solution to system of linear equations ---
inline void gesvx(char* fact, char* trans, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void gesvx(char* fact, char* trans, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void gesvx(char* fact, char* trans, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void gesvx(char* fact, char* trans, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equation ---
inline void gesvxx(char* fact, char* trans, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  sgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gesvxx(char* fact, char* trans, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void gesvxx(char* fact, char* trans, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  cgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void gesvxx(char* fact, char* trans, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes an LU factorization ---
inline void getrf2(integer* m, integer* n, float* a, integer* lda, integer* ipiv, integer* info)
{
  sgetrf2_(m, n, a, lda, ipiv, info);
}
inline void getrf2(integer* m, integer* n, double* a, integer* lda, integer* ipiv, integer* info)
{
  dgetrf2_(m, n, a, lda, ipiv, info);
}
inline void getrf2(integer* m, integer* n, scomplex* a, integer* lda, integer* ipiv, integer* info)
{
  cgetrf2_(m, n, a, lda, ipiv, info);
}
inline void getrf2(integer* m, integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* info)
{
  zgetrf2_(m, n, a, lda, ipiv, info);
}

// --- computes the inverse of a matrix using the LU factorization ---
inline void getri(integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* lwork, integer* info)
{
  sgetri_(n, a, lda, ipiv, work, lwork, info);
}
inline void getri(integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* lwork, integer* info)
{
  dgetri_(n, a, lda, ipiv, work, lwork, info);
}
inline void getri(integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  cgetri_(n, a, lda, ipiv, work, lwork, info);
}
inline void getri(integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zgetri_(n, a, lda, ipiv, work, lwork, info);
}

// --- solves a system of linear equations ---
inline void getrs(char* trans, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  sgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void getrs(char* trans, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  dgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void getrs(char* trans, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  cgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void getrs(char* trans, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves overdetermined or underdetermined real linear systems ---
inline void getsls(char* trans, integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  sgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline void getsls(char* trans, integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline void getsls(char* trans, integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  cgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline void getsls(char* trans, integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}

// --- forms the right or left eigenvectors of a real generalized eigenvalue problem ---
inline void ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* lscale,  float* rscale, integer* m, float* v, integer* ldv, integer* info)
{
  sggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}
inline void ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* lscale,  double* rscale, integer* m, double* v, integer* ldv, integer* info)
{
  dggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}
inline void ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* lscale,  float* rscale, integer* m, scomplex* v, integer* ldv, integer* info)
{
  cggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}
inline void ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* lscale,  double* rscale, integer* m, dcomplex* v, integer* ldv, integer* info)
{
  zggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}

// --- balances a pair of general real matrices (A,B) ---
inline void ggbal(char* job, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* ilo, integer* ihi, float* lscale, float* rscale, float* work, integer* info)
{
  sggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}
inline void ggbal(char* job, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* ilo, integer* ihi, double* lscale, double* rscale, double* work, integer* info)
{
  dggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}
inline void ggbal(char* job, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* ilo, integer* ihi, float* lscale, float* rscale, float* work, integer* info)
{
  cggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}
inline void ggbal(char* job, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* ilo, integer* ihi, double* lscale, double* rscale, double* work, integer* info)
{
  zggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}

// --- GGES computes the eigenvalues---
inline void gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* sdim, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* work, integer* lwork, logical* bwork, integer* info)
{
  sgges_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline void gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* sdim, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* work, integer* lwork, logical* bwork, integer* info)
{
  dgges_(jobvsl, jobvsr, sort, (L_fp) selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline void gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  cgges_(jobvsl, jobvsr, sort, (L_fp) selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}
inline void gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  zgges_(jobvsl, jobvsr, sort, (L_fp) selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}

inline void gegs(char* jobvsl, char* jobvsr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* work, integer* lwork, integer* info)
{
  printf(" Function sgegs() has been deprecated. Please use sgges() instead.\n");  
  sgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, info);
}
inline void gegs(char* jobvsl, char* jobvsr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* work, integer* lwork, integer* info)
{
  printf(" Function dgegs() has been deprecated. Please use dgges() instead.\n");  
  dgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, info);
}
inline void gegs(char* jobvsl, char* jobvsr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  printf(" Function cgegs() has been deprecated. Please use cgges() instead.\n");  
  cgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, info);
}
inline void gegs(char* jobvsl, char* jobvsr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  printf(" Function zgegs() has been deprecated. Please use zgges() instead.\n"); 
  zgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, info);
}

// --- SGGES3 computes the eigenvalues---
inline void gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* sdim, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* work, integer* lwork, logical* bwork, integer* info)
{
  sgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline void gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* sdim, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* work, integer* lwork, logical* bwork, integer* info)
{
  dgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline void gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  cgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}
inline void gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  zgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}

// --- computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices ---
inline void ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* sdim, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  sggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline void ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* sdim, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  dggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline void ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  cggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, info);
}
inline void ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  zggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, info);
}

// --- computes the eigenvalues, the left and/or right eigenvectors for GE matrices---
inline void ggev(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  sggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void ggev(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  dggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void ggev(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline void ggev(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

inline void gegv(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  printf(" Function sgegv() has been deprecated. Please use sggev() instead.\n"); 
  sgegv_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void gegv(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  printf(" Function dgegv() has been deprecated. Please use dggev() instead.\n");  
  dgegv_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void gegv(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  printf(" Function cgegv() has been deprecated. Please use cggev() instead.\n");  
  cgegv_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline void gegv(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  printf(" Function zgegv() has been deprecated. Please use zggev() instead.\n");  
  zgegv_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

// --- computes the eigenvalues, the left and/or right eigenvectors for GE matrices (blocked algorithm)---
inline void ggev3(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  sggev3_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void ggev3(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  dggev3_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void ggev3(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cggev3_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline void ggev3(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zggev3_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

// --- SGGEVX computes the eigenvalues, the left and/or right eigenvectors for GE matrices ---
inline void ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* ilo, integer* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, logical* bwork, integer* info)
{
  sggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info);
}
inline void ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* ilo, integer* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, logical* bwork, integer* info)
{
  dggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info);
}
inline void ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* ilo, integer* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, integer* iwork, logical* bwork, integer* info)
{
  cggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info);
}
inline void ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* ilo, integer* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, integer* iwork, logical* bwork, integer* info)
{
  zggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info);
}

// --- solves a general Gauss-Markov linear model (GLM) problem ---
inline void ggglm(integer* n, integer* m, integer* p, float* a, integer* lda, float* b, integer* ldb, float* d, float* x, float* y, float* work, integer* lwork, integer* info)
{
  sggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}
inline void ggglm(integer* n, integer* m, integer* p, double* a, integer* lda, double* b, integer* ldb, double* d, double* x, double* y, double* work, integer* lwork, integer* info)
{
  dggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}
inline void ggglm(integer* n, integer* m, integer* p, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* d, scomplex* x, scomplex* y, scomplex* work, integer* lwork, integer* info)
{
  cggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}
inline void ggglm(integer* n, integer* m, integer* p, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* d, dcomplex* x, dcomplex* y, dcomplex* work, integer* lwork, integer* info)
{
  zggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}

// --- reduces a pair of real matrices (A,B) ---
inline void gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* b, integer* ldb, float* q, integer* ldq, float* z, integer* ldz, float* work, integer* lwork, integer* info)
{
  sgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}
inline void gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* b, integer* ldb, double* q, integer* ldq, double* z, integer* ldz, double* work, integer* lwork, integer* info)
{
  dgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}
inline void gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* q, integer* ldq, scomplex* z, integer* ldz, scomplex* work, integer* lwork, integer* info)
{
  cgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}
inline void gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, integer* info)
{
  zgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}

// --- reduces a pair of real matrices (A,B) ---
inline void gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* b, integer* ldb, float* q, integer* ldq, float* z, integer* ldz, integer* info)
{
  sgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}
inline void gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* b, integer* ldb, double* q, integer* ldq, double* z, integer* ldz, integer* info)
{
  dgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}
inline void gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* q, integer* ldq, scomplex* z, integer* ldz, integer* info)
{
  cgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}
inline void gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, integer* info)
{
  zgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}

// --- solves overdetermined or underdetermined systems for OTHER matrices ---
inline void gglse(integer* m, integer* n, integer* p, float* a, integer* lda, float* b, integer* ldb, float* c, float* d, float* x, float* work, integer* lwork, integer* info)
{
  sgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}
inline void gglse(integer* m, integer* n, integer* p, double* a, integer* lda, double* b, integer* ldb, double* c, double* d, double* x, double* work, integer* lwork, integer* info)
{
  dgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}
inline void gglse(integer* m, integer* n, integer* p, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* c, scomplex* d, scomplex* x, scomplex* work, integer* lwork, integer* info)
{
  cgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}
inline void gglse(integer* m, integer* n, integer* p, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* c, dcomplex* d, dcomplex* x, dcomplex* work, integer* lwork, integer* info)
{
  zgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}

// --- computes a generalized QR factorization ---
inline void ggqrf(integer* n, integer* m, integer* p, float* a, integer* lda, float* taua, float* b, integer* ldb, float* taub, float* work, integer* lwork, integer* info)
{
  sggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline void ggqrf(integer* n, integer* m, integer* p, double* a, integer* lda, double* taua, double* b, integer* ldb, double* taub, double* work, integer* lwork, integer* info)
{
  dggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline void ggqrf(integer* n, integer* m, integer* p, scomplex* a, integer* lda, scomplex* taua, scomplex* b, integer* ldb, scomplex* taub, scomplex* work, integer* lwork, integer* info)
{
  cggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline void ggqrf(integer* n, integer* m, integer* p, dcomplex* a, integer* lda, dcomplex* taua, dcomplex* b, integer* ldb, dcomplex* taub, dcomplex* work, integer* lwork, integer* info)
{
  zggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}

// --- computes a generalized RQ factorization ---
inline void ggrqf(integer* m, integer* p, integer* n, float* a, integer* lda, float* taua, float* b, integer* ldb, float* taub, float* work, integer* lwork, integer* info)
{
  sggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline void ggrqf(integer* m, integer* p, integer* n, double* a, integer* lda, double* taua, double* b, integer* ldb, double* taub, double* work, integer* lwork, integer* info)
{
  dggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline void ggrqf(integer* m, integer* p, integer* n, scomplex* a, integer* lda, scomplex* taua, scomplex* b, integer* ldb, scomplex* taub, scomplex* work, integer* lwork, integer* info)
{
  cggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline void ggrqf(integer* m, integer* p, integer* n, dcomplex* a, integer* lda, dcomplex* taua, dcomplex* b, integer* ldb, dcomplex* taub, dcomplex* work, integer* lwork, integer* info)
{
  zggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}

// --- computes the singular value decomposition (SVD) for OTHER matrices ---
inline void ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, float* a, integer* lda, float* b, integer* ldb, float* alpha, float* beta, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, float* work, integer* lwork, integer* iwork, integer* info)
{
  sggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
}
inline void ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, double* a, integer* lda, double* b, integer* ldb, double* alpha, double* beta, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, double* work, integer* lwork, integer* iwork, integer* info)
{
  dggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
}
inline void ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* alpha, float* beta, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  cggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info);
}
inline void ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* alpha, double* beta, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  zggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info);
}

// --- computes the singular value decomposition ---
inline void ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, float* a, integer* lda, float* b, integer* ldb, float* alpha, float* beta, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, float* work, integer* iwork, integer* info)
{
  printf(" Function sggsvd() has been deprecated. Please use sggsvd3() instead.\n"); 
  sggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, iwork, info);
}
inline void ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, double* a, integer* lda, double* b, integer* ldb, double* alpha, double* beta, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, double* work, integer* iwork, integer* info)
{
  printf(" Function dggsvd() has been deprecated. Please use dggsvd3() instead.\n");  
  dggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, iwork, info);
}
inline void ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* alpha, float* beta, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, scomplex* work, float* rwork, integer* iwork, integer* info)
{
  printf(" Function cggsvd() has been deprecated. Please use cggsvd3() instead.\n");  
  cggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, rwork, iwork, info);
}
inline void ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* alpha, double* beta, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, dcomplex* work, double* rwork, integer* iwork, integer* info)
{
  printf(" Function dggsvd() has been deprecated. Please use dggsvd3() instead.\n");  
  zggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, rwork, iwork, info);
}

// --- computes orthogonal matrices U, V and Q ---
inline void ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, float* a, integer* lda, float* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, integer* iwork, float* tau, float* work, integer* lwork, integer* info)
{
  sggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info);
}
inline void ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, double* a, integer* lda, double* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, integer* iwork, double* tau, double* work, integer* lwork, integer* info)
{
  dggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info);
}
inline void ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, integer* iwork, float* rwork, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info);
}
inline void ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, integer* iwork, double* rwork, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info);
}

// --- computes orthogonal matrices U, V and Q ---
inline void ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, float* a, integer* lda, float* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, integer* iwork, float* tau, float* work, integer* info)
{
  printf(" Function sggsvp() has been deprecated. Please use sggsvp3() instead.\n");  
  sggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, info);
}
inline void ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, double* a, integer* lda, double* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, integer* iwork, double* tau, double* work, integer* info)
{
  printf(" Function dggsvp() has been deprecated. Please use dggsvp3() instead.\n");  
  dggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, info);
}
inline void ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, integer* iwork, float* rwork, scomplex* tau, scomplex* work, integer* info)
{
  printf(" Function cggsvp() has been deprecated. Please use cggsvp3() instead.\n");  
  cggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, info);
}
inline void ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, integer* iwork, double* rwork, dcomplex* tau, dcomplex* work, integer* info)
{
  printf(" Function zggsvp() has been deprecated. Please use zggsvp3() instead.\n"); 
  zggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, info);
}

// --- estimates the reciprocal of the condition number ---
inline void gtcon(char* norm, integer* n,  float* dl,  float* d, float* du, float* du2, integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  sgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, iwork, info);
}
inline void gtcon(char* norm, integer* n,  double* dl, double* d, double* du, double* du2, integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, iwork, info);
}
inline void gtcon(char* norm, integer* n,  scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  cgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, info);
}
inline void gtcon(char* norm, integer* n,  dcomplex* dl, dcomplex* d,  dcomplex* du, dcomplex* du2, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  zgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, info);
}

// --- improves the computed solution to a system of linear equations ---
inline void gtrfs(char* trans, integer* n, integer* nrhs, float* dl, float* d,  float* du,  float* dlf, float* df, float* duf, float* du2, integer* ipiv, float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void gtrfs(char* trans, integer* n, integer* nrhs,  double* dl,  double* d,  double* du,  double* dlf,  double* df,  double* duf,  double* du2,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void gtrfs(char* trans, integer* n, integer* nrhs,  scomplex* dl,  scomplex* d,  scomplex* du, scomplex* dlf,  scomplex* df,  scomplex* duf,  scomplex* du2,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void gtrfs(char* trans, integer* n, integer* nrhs,  dcomplex* dl,  dcomplex* d,  dcomplex* du,  dcomplex* dlf,  dcomplex* df,  dcomplex* duf,  dcomplex* du2,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations ---
inline void gtsv(integer* n, integer* nrhs, float* dl, float* d, float* du, float* b, integer* ldb, integer* info)
{
  sgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}
inline void gtsv(integer* n, integer* nrhs, double* dl, double* d, double* du, double* b, integer* ldb, integer* info)
{
  dgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}
inline void gtsv(integer* n, integer* nrhs, scomplex* dl, scomplex* d, scomplex* du, scomplex* b, integer* ldb, integer* info)
{
  cgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}
inline void gtsv(integer* n, integer* nrhs, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* b, integer* ldb, integer* info)
{
  zgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}

// --- uses the LU factorization to compute the solution ---
inline void gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  float* dl,  float* d, float* du, float* dlf, float* df, float* duf, float* du2, integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  double* dl,  double* d,  double* du, double* dlf, double* df, double* duf, double* du2, integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  scomplex* dl,  scomplex* d,  scomplex* du, scomplex* dlf, scomplex* df, scomplex* duf, scomplex* du2, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  dcomplex* dl,  dcomplex* d,  dcomplex* du, dcomplex* dlf, dcomplex* df, dcomplex* duf, dcomplex* du2, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes an LU factorization of a real tridiagonal matrix A ---
inline void gttrf(integer* n, float* dl, float* d, float* du, float* du2, integer* ipiv, integer* info)
{
  sgttrf_(n, dl, d, du, du2, ipiv, info);
}
inline void gttrf(integer* n, double* dl, double* d, double* du, double* du2, integer* ipiv, integer* info)
{
  dgttrf_(n, dl, d, du, du2, ipiv, info);
}
inline void gttrf(integer* n, scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, integer* ipiv, integer* info)
{
  cgttrf_(n, dl, d, du, du2, ipiv, info);
}
inline void gttrf(integer* n, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* du2, integer* ipiv, integer* info)
{
  zgttrf_(n, dl, d, du, du2, ipiv, info);
}

// --- solves one of the systems of equations ---
inline void gttrs(char* trans, integer* n, integer* nrhs,  float* dl,  float* d,  float* du,  float* du2,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  sgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}
inline void gttrs(char* trans, integer* n, integer* nrhs,  double* dl,  double* d,  double* du,  double* du2,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  dgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}
inline void gttrs(char* trans, integer* n, integer* nrhs,  scomplex* dl,  scomplex* d,  scomplex* du,  scomplex* du2,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  cgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}
inline void gttrs(char* trans, integer* n, integer* nrhs,  dcomplex* dl,  dcomplex* d,  dcomplex* du,  dcomplex* du2,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}

// --- estimates the 1-norm of a square matrix ---
inline void lacn2(integer* n, float* v, float* x, integer* isgn, float* est, integer* kase, integer* isave)
{
  slacn2_(n, v, x, isgn, est, kase, isave);
}
inline void lacn2(integer* n, double* v, double* x, integer* isgn, double* est, integer* kase, integer* isave)
{
  dlacn2_(n, v, x, isgn, est, kase, isave);
}
inline void lacn2(integer* n, scomplex* v, scomplex* x, float* est, integer* kase, integer* isave)
{
  clacn2_(n, v, x, est, kase, isave);
}
inline void lacn2(integer* n, dcomplex* v, dcomplex* x, double* est, integer* kase, integer* isave)
{
  zlacn2_(n, v, x, est, kase, isave);
}

// --- copies all or part of one two-dimensional array ---
inline void lacpy(char* uplo, integer* m, integer* n,  float* a, integer* lda, float* b, integer* ldb)
{
  slacpy_(uplo, m, n,  a, lda, b, ldb);
}
inline void lacpy(char* uplo, integer* m, integer* n,  double* a, integer* lda, double* b, integer* ldb)
{
  dlacpy_(uplo, m, n,  a, lda, b, ldb);
}
inline void lacpy(char* uplo, integer* m, integer* n,  scomplex* a, integer* lda, scomplex* b, integer* ldb)
{
  clacpy_(uplo, m, n,  a, lda, b, ldb);
}
inline void lacpy(char* uplo, integer* m, integer* n, dcomplex* a, integer* lda,dcomplex* b, integer* ldb)
{
  zlacpy_(uplo, m, n,  a, lda, b, ldb);
}

inline void slag2d(integer* m, integer* n, float* sa, integer* ldsa, double* a, integer* lda, integer* info)
{
  slag2d_(m, n, sa, ldsa, a, lda, info);
}

// --- converts a double precision matrix to a single precision matrix ---
inline void dlag2s(integer*m, integer* n, double* a, integer*lda, float* sa, integer* ldsa, integer* info)
{
  dlag2s_(m, n, a, lda, sa, ldsa, info);
}
// --- converts a complex single precision matrix to a complex double precision matrix. ---
inline void clag2z(integer*m, integer* n, scomplex* sa, integer*ldsa, dcomplex* a, integer* lda, integer* info)
{
  clag2z_(m, n, sa, ldsa, a, lda, info);
} 
// --- converts a single precision matrix to a double precision matrix ---
inline void zlag2c(integer*m, integer* n, dcomplex* a, integer*lda, scomplex* sa, integer* ldsa, integer* info)
{
  zlag2c_(m, n, a, lda, sa, ldsa, info);
}

// --- returns sqrt(x2+y2) ---
inline float lapy2(float* x, float* y)
{
  return slapy2_(x, y);
}
inline double lapy2(double* x, double* y)
{
  return dlapy2_(x, y);
}

// --- returns sqrt(x2+y2+z2) ---
inline float lapy3(float* x, float* y, float* z)
{
  return slapy3_(x, y, z);
}
inline double lapy3(double* x, double* y, double* z)
{
  return dlapy3_(x, y, z);
}

// --- generates a plane rotation so that the diagonal is nonnegative ---
inline void lartgp(float* f, float* g, float* cs, float* sn, float* r)
{
  slartgp_(f, g, cs, sn, r);
}
inline void lartgp(double* f, double* g, double* cs, double* sn, double* r)
{
  dlartgp_(f, g, cs, sn, r);
}

// --- generates a plane rotation designed to introduce a bulge in implicit QR iteration for the bidiagonal SVD problem ---
inline void lartgs(float* x, float* y, float* sigma, float* cs, float* sn)
{
  slartgs_(x, y, sigma, cs, sn);
}
inline void lartgs(double* x, double* y, double* sigma, double* cs, double* sn)
{
  dlartgs_(x, y, sigma, cs, sn);
}

// ---  returns the value of the largest absolute value of any element of a general rectangular matrix ---
inline float lange(char* norm, integer* m, integer* n,  float* a, integer* lda, float* work)
{
  return slange_(norm, m, n,  a, lda, work);
}
inline double lange(char* norm, integer* m, integer* n,  double* a, integer* lda, double* work)
{
  return dlange_(norm, m, n,  a, lda, work);
}
inline float lange(char* norm, integer* m, integer* n,  scomplex* a, integer* lda, float* work)
{
  return clange_(norm, m, n,  a, lda, work);
}
inline double lange(char* norm, integer* m, integer* n,  dcomplex* a, integer* lda, double* work)
{
  return zlange_(norm, m, n,  a, lda, work);
}

// --- returns the largest absolute value of a real symmetric matrix---
inline float lansy(char* norm, char* uplo, integer* n,  float* a, integer* lda, float* work)
{
  return slansy_(norm, uplo, n,  a, lda, work);
}
inline double lansy(char* norm, char* uplo, integer* n,  double* a, integer* lda, double* work)
{
  return dlansy_(norm, uplo, n,  a, lda, work);
}
inline float lansy(char* norm, char* uplo, integer* n,  scomplex* a, integer* lda, float* work)
{
  return clansy_(norm, uplo, n,  a, lda, work);
}
inline double lansy(char* norm, char* uplo, integer* n,  dcomplex* a, integer* lda, double* work)
{
  return zlansy_(norm, uplo, n,  a, lda, work);
}

// --- returns the largest absolute value of a trapezoidal or triangular matrix---
inline float lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  float* a, integer* lda, float* work)
{
  return slantr_(norm, uplo, diag, m, n,  a, lda, work);
}
inline double lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  double* a, integer* lda, double* work)
{
  return dlantr_(norm, uplo, diag, m, n,  a, lda, work);
}
inline float lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  scomplex* a, integer* lda, float* work)
{
  return clantr_(norm, uplo, diag, m, n,  a, lda, work);
}
inline double lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  dcomplex* a, integer* lda, double* work)
{
  return zlantr_(norm, uplo, diag, m, n,  a, lda, work);
}

// --- rearranges rows of a matrix as specified by a permutation vector. ---
inline void lapmr(logical* forwrd, integer* m, integer* n, float* x, integer* ldx, integer* k)
{
  slapmr_(forwrd, m, n, x, ldx, k);
}
inline void lapmr(logical* forwrd, integer* m, integer* n, double* x, integer* ldx, integer* k)
{
  dlapmr_(forwrd, m, n, x, ldx, k);
}
inline void lapmr(logical* forwrd, integer* m, integer* n, scomplex* x, integer* ldx, integer* k)
{
  clapmr_(forwrd, m, n, x, ldx, k);
}
inline void lapmr(logical* forwrd, integer* m, integer* n, dcomplex* x, integer* ldx, integer* k)
{
  zlapmr_(forwrd, m, n, x, ldx, k);
}

// --- performs a forward or backward permutation of the columns of a matrix ---
inline void lapmt(logical* forwrd, integer* m, integer* n, float* x, integer* ldx, integer* k)
{
  slapmt_(forwrd, m, n, x, ldx, k);
}
inline void lapmt(logical* forwrd, integer* m, integer* n, double* x, integer* ldx, integer* k)
{
  dlapmt_(forwrd, m, n, x, ldx, k);
}
inline void lapmt(logical* forwrd, integer* m, integer* n, scomplex* x, integer* ldx, integer* k)
{
  clapmt_(forwrd, m, n, x, ldx, k);
}
inline void lapmt(logical* forwrd, integer* m, integer* n, dcomplex* x, integer* ldx, integer* k)
{
  zlapmt_(forwrd, m, n, x, ldx, k);
}

// --- applies a block reflector or its transpose to a general rectangular matrix ---
inline void larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  float* v, integer* ldv,  float* t, integer* ldt, float* c, integer* ldc, float* work, integer* ldwork)
{
  slarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline void larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  double* v, integer* ldv,  double* t, integer* ldt, double* c, integer* ldc, double* work, integer* ldwork)
{
  dlarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline void larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* c, integer* ldc, scomplex* work, integer* ldwork)
{
  clarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline void larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* c, integer* ldc, dcomplex* work, integer* ldwork)
{
  zlarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}

// --- applies an elementary reflector to a general rectangular matrix ---
inline void larfx(char* side, integer* m, integer* n, float* v, float* tau, float* c, integer* ldc, float* work)
{
  slarfx_(side, m, n,  v, tau, c, ldc, work);
}
inline void larfx(char* side, integer* m, integer* n, double* v, double* tau, double* c, integer* ldc, double* work)
{
  dlarfx_(side, m, n,  v, tau, c, ldc, work);
}
inline void larfx(char* side, integer* m, integer* n, scomplex* v, scomplex* tau, scomplex* c, integer* ldc, scomplex* work)
{
  clarfx_(side, m, n,  v, tau, c, ldc, work);
}
inline void larfx(char* side, integer* m, integer* n, dcomplex* v, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work)
{
  zlarfx_(side, m, n,  v, tau, c, ldc, work);
}

// ---  returns a vector of random numbers from a uniform or normal distribution ---
inline void larnv(integer* idist, integer* iseed, integer* n, float* x)
{
  slarnv_(idist, iseed, n, x);
}
inline void larnv(integer* idist, integer* iseed, integer* n, double* x)
{
  dlarnv_(idist, iseed, n, x);
}
inline void larnv(integer* idist, integer* iseed, integer* n, scomplex* x)
{
  clarnv_(idist, iseed, n, x);
}
inline void larnv(integer* idist, integer* iseed, integer* n, dcomplex* x)
{
  zlarnv_(idist, iseed, n, x);
}

// --- multiplies a general rectangular matrix by a real scalar defined as cto/cfrom ---
inline void lascl(char* type, integer* kl, integer* ku, float* cfrom, float* cto, integer* m, integer* n, float* a, integer* lda, integer* info)
{
  slascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}
inline void lascl(char* type, integer* kl, integer* ku, double* cfrom, double* cto, integer* m, integer* n, double* a, integer* lda, integer* info)
{
  dlascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}
inline void lascl(char* type, integer* kl, integer* ku, float* cfrom, float* cto, integer* m, integer* n, scomplex* a, integer* lda, integer* info)
{
  clascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}
inline void lascl(char* type, integer* kl, integer* ku, double* cfrom, double* cto, integer* m, integer* n, dcomplex* a, integer* lda, integer* info)
{
  zlascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}

// --- sorts numbers in increasing or decreasing order. ---
inline void lasrt(char* id, integer* n, float* d, integer* info)
{
  slasrt_(id, n, d, info);
}
inline void lasrt(char* id, integer* n, double* d, integer* info)
{
  dlasrt_(id, n, d, info);
}

// --- updates a sum of squares represented in scaled form ---
inline void lassq(integer* n, float* x, integer* incx, float* scale, float* sumsq)
{
  slassq_(n, x, incx, scale, sumsq);
}
inline void lassq(integer* n, double* x, integer* incx, double* scale, double* sumsq)
{
  dlassq_(n, x, incx, scale, sumsq);
}
inline void lassq(integer* n, scomplex* x, integer* incx, float* scale, float* sumsq)
{
  classq_(n, x, incx, scale, sumsq);
}
inline void lassq(integer* n, dcomplex* x, integer* incx, double* scale, double* sumsq)
{
  zlassq_(n, x, incx, scale, sumsq);
}

// --- generates a real orthogonal matrix Q ---
inline void opgtr(char* uplo, integer* n, float* ap, float* tau, float* q, integer* ldq, float *work, integer *info)
{
  sopgtr_(uplo, n, ap, tau, q, ldq, work, info);
}
inline void opgtr(char* uplo, integer* n, double* ap, double* tau, double* q, integer* ldq, double *work, integer *info)
{
  dopgtr_(uplo, n, ap, tau, q, ldq, work, info);
}
inline void upgtr(char* uplo, integer* n, scomplex* ap, scomplex* tau, scomplex* q, integer* ldq, scomplex* work, integer* info)
{
  cupgtr_(uplo, n, ap, tau, q, ldq, work, info);
}
inline void upgtr(char* uplo, integer* n, dcomplex* ap, dcomplex* tau, dcomplex* q, integer* ldq, dcomplex* work, integer* info)
{
  zupgtr_(uplo, n, ap, tau, q, ldq, work, info);
}

// --- overwrites the general real M-by-N matrix ---
inline void opmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  float* ap,  float* tau, float* c, integer* ldc, float* work, integer* info)
{
  sopmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}
inline void opmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  double* ap,  double* tau, double* c, integer* ldc, double* work, integer* info)
{
  dopmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}

inline void upmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  scomplex* ap, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  cupmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}
inline void upmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  dcomplex* ap,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  zupmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}

// --- simultaneously bidiagonalizes the blocks of an M-by-M partitioned orthogonal matrix X ---
inline void orbdb(char* trans, char* signs, integer* m, integer* p, integer* q, float* x11, integer* ldx11, float* x12, integer* ldx12, float* x21, integer* ldx21, float* x22, integer* ldx22, float* theta, float* phi, float* taup1, float* taup2, float* tauq1, float* tauq2, float* work, integer* lwork, integer* info)
{
  sorbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}
inline void orbdb(char* trans, char* signs, integer* m, integer* p, integer* q, double* x11, integer* ldx11, double* x12, integer* ldx12, double* x21, integer* ldx21, double* x22, integer* ldx22, double* theta, double* phi, double* taup1, double* taup2, double* tauq1, double* tauq2, double* work, integer* lwork, integer* info)
{
  dorbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}
inline void unbdb(char* trans, char* signs, integer* m, integer* p, integer* q, scomplex* x11, integer* ldx11, scomplex* x12, integer* ldx12, scomplex* x21, integer* ldx21, scomplex* x22, integer* ldx22, float* theta, float* phi, scomplex* taup1, scomplex* taup2, scomplex* tauq1, scomplex* tauq2, scomplex* work, integer* lwork, integer* info)
{
  cunbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}
inline void unbdb(char* trans, char* signs, integer* m, integer* p, integer* q, dcomplex* x11, integer* ldx11, dcomplex* x12, integer* ldx12, dcomplex* x21, integer* ldx21, dcomplex* x22, integer* ldx22, double* theta, double* phi, dcomplex* taup1, dcomplex* taup2, dcomplex* tauq1, dcomplex* tauq2, dcomplex* work, integer* lwork, integer* info)
{
  zunbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}

// --- computes the CS decomposition of an M-by-M partitioned orthogonal matrix X ---
inline void orcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, float* x11, integer* ldx11, float* x12, integer* ldx12, float* x21, integer* ldx21, float* x22, integer* ldx22, float* theta, float* u1, integer* ldu1, float* u2, integer* ldu2, float* v1t, integer* ldv1t, float* v2t, integer* ldv2t, float *work, integer *lwork, integer *iwork, integer *info)
{
  sorcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, iwork, info);
}
inline void orcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, double* x11, integer* ldx11, double* x12, integer* ldx12, double* x21, integer* ldx21, double* x22, integer* ldx22, double* theta, double* u1, integer* ldu1, double* u2, integer* ldu2, double* v1t, integer* ldv1t, double* v2t, integer* ldv2t, double *work, integer *lwork, integer *iwork, integer *info)
{
  dorcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, iwork, info);
}
inline void uncsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, scomplex* x11, integer* ldx11, scomplex* x12, integer* ldx12, scomplex* x21, integer* ldx21, scomplex* x22, integer* ldx22, float* theta, scomplex* u1, integer* ldu1, scomplex* u2, integer* ldu2, scomplex* v1t, integer* ldv1t, scomplex* v2t, integer* ldv2t, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* info)
{
  cuncsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, rwork, lrwork, iwork, info);
}
inline void uncsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, dcomplex* x11, integer* ldx11, dcomplex* x12, integer* ldx12, dcomplex* x21, integer* ldx21, dcomplex* x22, integer* ldx22, double* theta, dcomplex* u1, integer* ldu1, dcomplex* u2, integer* ldu2, dcomplex* v1t, integer* ldv1t, dcomplex* v2t, integer* ldv2t, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* info)
{
  zuncsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, rwork, lrwork, iwork, info);
}

// --- computes the CS decomposition of an M-by-Q matrix X with orthonormal columns ---
inline void orcsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, float* x11, integer* ldx11, float* x21, integer* ldx21, float* theta, float* u1, integer* ldu1, float* u2, integer* ldu2, float* v1t, integer* ldv1t, float *work, integer *lwork, integer *iwork, integer *info)
{
  sorcsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, iwork, info);
}
inline void orcsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, double* x11, integer* ldx11, double* x21, integer* ldx21, double* theta, double* u1, integer* ldu1, double* u2, integer* ldu2, double* v1t, integer* ldv1t, double *work, integer *lwork, integer *iwork, integer *info)
{
  dorcsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, iwork, info);
}
inline void uncsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, scomplex* x11, integer* ldx11, scomplex* x21, integer* ldx21, float* theta, scomplex* u1, integer* ldu1, scomplex* u2, integer* ldu2, scomplex* v1t, integer* ldv1t, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* info)
{
  cuncsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, rwork, lrwork, iwork, info);
}
inline void uncsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, dcomplex* x11, integer* ldx11, dcomplex* x21, integer* ldx21, double* theta, dcomplex* u1, integer* ldu1, dcomplex* u2, integer* ldu2, dcomplex* v1t, integer* ldv1t, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* info)
{
  zuncsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, rwork, lrwork, iwork, info);
}

// --- generates an M-by-N real matrix Q with orthonormal columns ---
inline void orgql(integer* m, integer* n, integer* k, float* a, integer* lda,  float* tau, float* work, integer* lwork, integer* info)
{
  sorgql_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void orgql(integer* m, integer* n, integer* k, double* a, integer* lda,  double* tau, double* work, integer* lwork, integer* info)
{
  dorgql_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void ungql(integer* m, integer* n, integer* k, scomplex* a, integer* lda,  scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cungql_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void ungql(integer* m, integer* n, integer* k, dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zungql_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- generates an M-by-N real matrix Q with orthonormal rows ---
inline void orgrq(integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  sorgrq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void orgrq(integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dorgrq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void ungrq(integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cungrq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline void ungrq(integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zungrq_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- overwrites the general real M-by-N matrix ---
inline void ormql(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sormql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline void ormql(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dormql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline void unmql(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline void unmql(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}

// --- overwrites the general real M-by-N matrix ---
inline void ormrq(char* side, char* trans, integer* m, integer* n, integer* k,  float* a, integer* lda,  float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sormrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline void ormrq(char* side, char* trans, integer* m, integer* n, integer* k,  double* a, integer* lda,  double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dormrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmrq(char* side, char* trans, integer* m, integer* n, integer* k,  scomplex* a, integer* lda,  scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmrq(char* side, char* trans, integer* m, integer* n, integer* k,  dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}

// --- overwrites the general real M-by-N matrix---
inline void ormrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  float* a, integer* lda,  float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  sormrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}
inline void ormrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  double* a, integer* lda,  double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  dormrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  scomplex* a, integer* lda,  scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}
inline void unmrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}

inline void latzm(char* side, integer* m, integer* n, float* v, integer* incv, float* tau, float* c1, float* c2, integer* ldc, float* work)
{
  slatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}
inline void latzm(char* side, integer* m, integer* n, double* v, integer* incv, double* tau, double* c1, double* c2, integer* ldc, double* work)
{
  dlatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}
inline void latzm(char* side, integer* m, integer* n, scomplex* v, integer* incv, scomplex* tau, scomplex* c1, scomplex* c2, integer* ldc, scomplex* work)
{
  printf(" Function clatzm() has been deprecated. Please use cunmrz() instead.\n"); 
  clatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}
inline void latzm(char* side, integer* m, integer* n, dcomplex* v, integer* incv, dcomplex* tau, dcomplex* c1, dcomplex* c2, integer* ldc, dcomplex* work)
{
  printf(" Function zlatzm() has been deprecated. Please use zunmrz() instead.\n"); 
  zlatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}

// --- generates a real orthogonal matrix Q ---
inline void orghr(integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* tau, float *work, integer *lwork, integer *info)
{
  sorghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}
inline void orghr(integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* tau, double *work, integer *lwork, integer *info)
{
  dorghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}
inline void unghr(integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  cunghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}
inline void unghr(integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zunghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}

// --- overwrites the general complex M-by-N matrix C ---
inline void ormhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  float* a, integer* lda,  float* tau, float* c, integer* ldc, float *work, integer *lwork, integer *info)
{
  sormhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline void ormhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  double* a, integer* lda,  double* tau, double* c, integer* ldc, double *work, integer *lwork, integer *info)
{
  dormhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline void unmhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  scomplex* a, integer* lda,  scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  cunmhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline void unmhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  zunmhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}

// --- estimates the reciprocal of the condition number ---
inline void pbcon(char* uplo, integer* n, integer* kd,  float* ab, integer* ldab, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  spbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, iwork, info);
}
inline void pbcon(char* uplo, integer* n, integer* kd,  double* ab, integer* ldab, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dpbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, iwork, info);
}
inline void pbcon(char* uplo, integer* n, integer* kd,  scomplex* ab, integer* ldab, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  cpbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, rwork, info);
}
inline void pbcon(char* uplo, integer* n, integer* kd,  dcomplex* ab, integer* ldab, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  zpbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings---
inline void pbequ(char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* s, float* scond, float* amax, integer* info)
{
  spbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}
inline void pbequ(char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* s, double* scond, double* amax, integer* info)
{
  dpbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}
inline void pbequ(char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* s, float* scond, float* amax, integer* info)
{
  cpbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}
inline void pbequ(char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* s, double* scond, double* amax, integer* info)
{
  zpbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}

// --- improves the computed solution to a system of linear equations ---
inline void pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab,  float* afb, integer* ldafb,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  spbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab,  double* afb, integer* ldafb,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dpbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab,  scomplex* afb, integer* ldafb,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cpbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab,  dcomplex* afb, integer* ldafb,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zpbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes a split Cholesky factorization of a real symmetric positive definite band matrix---
inline void pbstf(char* uplo, integer* n, integer* kb, float* bb, integer* ldbb, integer* info)
{
  spbstf_(uplo, n, kb, bb, ldbb, info);
}
inline void pbstf(char* uplo, integer* n, integer* kb, double* bb, integer* ldbb, integer* info)
{
  dpbstf_(uplo, n, kb, bb, ldbb, info);
}
inline void pbstf(char* uplo, integer* n, integer* kb, scomplex* bb, integer* ldbb, integer* info)
{
  cpbstf_(uplo, n, kb, bb, ldbb, info);
}
inline void pbstf(char* uplo, integer* n, integer* kb, dcomplex* bb, integer* ldbb, integer* info)
{
  zpbstf_(uplo, n, kb, bb, ldbb, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline void pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, float* ab, integer* ldab, float* b, integer* ldb, integer* info)
{
  spbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline void pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, double* ab, integer* ldab, double* b, integer* ldb, integer* info)
{
  dpbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline void pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, scomplex* ab, integer* ldab, scomplex* b, integer* ldb, integer* info)
{
  cpbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline void pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* b, integer* ldb, integer* info)
{
  zpbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}

// --- computes the solution to system of linear equations ---
inline void pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  spbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, double* ab, integer* ldab, double* afb, integer* ldafb, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dpbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cpbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zpbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the Cholesky factorization ---
inline void pbtrf(char* uplo, integer* n, integer* kd, float* ab, integer*ldab, integer* info)
{
  spbtrf_(uplo, n, kd, ab, ldab, info);
}
inline void pbtrf(char* uplo, integer* n, integer* kd, double* ab, integer*ldab, integer* info)
{
  dpbtrf_(uplo, n, kd, ab, ldab, info);
}
inline void pbtrf(char* uplo, integer* n, integer* kd, scomplex* ab, integer*ldab, integer* info)
{
  cpbtrf_(uplo, n, kd, ab, ldab, info);
}
inline void pbtrf(char* uplo, integer* n, integer* kd, dcomplex* ab, integer*ldab, integer* info)
{
  zpbtrf_(uplo, n, kd, ab, ldab, info);
}

// --- solves a system of linear equations ---
inline void pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab, float* b, integer* ldb, integer* info)
{
  spbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline void pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab, double* b, integer* ldb, integer* info)
{
  dpbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline void pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab, scomplex* b, integer* ldb, integer* info)
{
  cpbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline void pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab, dcomplex* b, integer* ldb, integer* info)
{
  zpbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}

// --- computes the Cholesky factorization of a real symmetric ---
inline void pftrf(char* transr, char* uplo, integer* n, float* a, integer* info)
{
  spftrf_(transr, uplo, n, a, info);
}
inline void pftrf(char* transr, char* uplo, integer* n, double* a, integer* info)
{
  dpftrf_(transr, uplo, n, a, info);
}
inline void pftrf(char* transr, char* uplo, integer* n, scomplex* a, integer* info)
{
  cpftrf_(transr, uplo, n, a, info);
}
inline void pftrf(char* transr, char* uplo, integer* n, dcomplex* a, integer* info)
{
  zpftrf_(transr, uplo, n, a, info);
}

// --- computes the inverse of a real (symmetric) positive definite matrix ---
inline void pftri(char* transr, char* uplo, integer* n, float* a, integer* info)
{
  spftri_(transr, uplo, n, a, info);
}
inline void pftri(char* transr, char* uplo, integer* n, double* a, integer* info)
{
  dpftri_(transr, uplo, n, a, info);
}
inline void pftri(char* transr, char* uplo, integer* n, scomplex* a, integer* info)
{
  cpftri_(transr, uplo, n, a, info);
}
inline void pftri(char* transr, char* uplo, integer* n, dcomplex* a, integer* info)
{
  zpftri_(transr, uplo, n, a, info);
}

// --- solves a system of linear equations A*X = B with a symmetric matrix ---
inline void pftrs(char* transr, char* uplo, integer* n, integer* nrhs,  float* a, float* b, integer* ldb, integer* info)
{
  spftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}
inline void pftrs(char* transr, char* uplo, integer* n, integer* nrhs,  double* a, double* b, integer* ldb, integer* info)
{
  dpftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}
inline void pftrs(char* transr, char* uplo, integer* n, integer* nrhs, scomplex* a, scomplex* b, integer* ldb, integer* info)
{
  cpftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}
inline void pftrs(char* transr, char* uplo, integer* n, integer* nrhs, dcomplex* a, dcomplex* b, integer* ldb, integer* info)
{
  zpftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}

// --- estimates the reciprocal of the condition number ---
inline void pocon(char* uplo, integer* n,  float* a, integer* lda, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  spocon_(uplo, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline void pocon(char* uplo, integer* n,  double* a, integer* lda, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dpocon_(uplo, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline void pocon(char* uplo, integer* n,  scomplex* a, integer* lda, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  cpocon_(uplo, n,  a, lda, anorm, rcond, work, rwork, info);
}
inline void pocon(char* uplo, integer* n,  dcomplex* a, integer* lda, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  zpocon_(uplo, n,  a, lda, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings ---
inline void poequ(integer* n,  float* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  spoequ_(n, a, lda, s, scond, amax, info);
}
inline void poequ(integer* n,  double* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  dpoequ_(n, a, lda, s, scond, amax, info);
}
inline void poequ(integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  cpoequ_(n, a, lda, s, scond, amax, info);
}
inline void poequ(integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  zpoequ_(n, a, lda, s, scond, amax, info);
}

// --- computes row and column scalings ---
inline void poequb(integer* n,  float* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  spoequb_(n, a, lda, s, scond, amax, info);
}
inline void poequb(integer* n,  double* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  dpoequb_(n, a, lda, s, scond, amax, info);
}
inline void poequb(integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  cpoequb_(n, a, lda, s, scond, amax, info);
}
inline void poequb(integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  zpoequb_(n, a, lda, s, scond, amax, info);
}

// --- improves the computed solution to a system of linear equations ---
inline void porfs(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void porfs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void porfs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void porfs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- improves the computed solution to a system of linear equations ---
inline void porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  float* s,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  sporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  double* s,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  float* s,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  cporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  double* s,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for PO matrices ---
inline void posv(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  sposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}
inline void posv(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  dposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}
inline void posv(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  cposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}
inline void posv(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  zposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for PO matrices ---
inline void posvx(char* fact, char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void posvx(char* fact, char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void posvx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void posvx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for PO matrices ---
inline void posvxx(char* fact, char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  sposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void posvxx(char* fact, char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void posvxx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  cposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void posvxx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// ---computes the Cholesky factorization of a real symmetric matrix A---
inline void potrf2(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  spotrf2_(uplo, n, a, lda, info);
}
inline void potrf2(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  dpotrf2_(uplo, n, a, lda, info);
}
inline void potrf2(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  cpotrf2_(uplo, n, a, lda, info);
}
inline void potrf2(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  zpotrf2_(uplo, n, a, lda, info);
}

// --- solves a system of linear equations A*X = B ---
inline void potrs(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  spotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}
inline void potrs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  dpotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}
inline void potrs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  cpotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}
inline void potrs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  zpotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}

// --- estimates the reciprocal of the condition number ---
inline void ppcon(char* uplo, integer* n, float* ap, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  sppcon_(uplo, n,  ap, anorm, rcond, work, iwork, info);
}
inline void ppcon(char* uplo, integer* n, double* ap, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dppcon_(uplo, n,  ap, anorm, rcond, work, iwork, info);
}
inline void ppcon(char* uplo, integer* n, scomplex* ap, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  cppcon_(uplo, n,  ap, anorm, rcond, work, rwork, info);
}
inline void ppcon(char* uplo, integer* n, dcomplex* ap, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  zppcon_(uplo, n,  ap, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings ---
inline void ppequ(char* uplo, integer* n, float* ap, float* s, float* scond, float* amax, integer* info)
{
  sppequ_(uplo, n, ap, s, scond, amax, info);
}
inline void ppequ(char* uplo, integer* n, double* ap, double* s, double* scond, double* amax, integer* info)
{
  dppequ_(uplo, n, ap, s, scond, amax, info);
}
inline void ppequ(char* uplo, integer* n, scomplex* ap, float* s, float* scond, float* amax, integer* info)
{
  cppequ_(uplo, n, ap, s, scond, amax, info);
}
inline void ppequ(char* uplo, integer* n, dcomplex* ap, double* s, double* scond, double* amax, integer* info)
{
  zppequ_(uplo, n, ap, s, scond, amax, info);
}

// --- improves the computed solution to a system of linear equations ---
inline void pprfs(char* uplo, integer* n, integer* nrhs,  float* ap,  float* afp,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  spprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void pprfs(char* uplo, integer* n, integer* nrhs,  double* ap,  double* afp,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dpprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void pprfs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  scomplex* afp,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cpprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void pprfs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* afp,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zpprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline void ppsv(char* uplo, integer* n, integer* nrhs, float* ap, float* b, integer* ldb, integer* info)
{
  sppsv_(uplo, n, nrhs, ap, b, ldb, info);
}
inline void ppsv(char* uplo, integer* n, integer* nrhs, double* ap, double* b, integer* ldb, integer* info)
{
  dppsv_(uplo, n, nrhs, ap, b, ldb, info);
}
inline void ppsv(char* uplo, integer* n, integer* nrhs, scomplex* ap, scomplex* b, integer* ldb, integer* info)
{
  cppsv_(uplo, n, nrhs, ap, b, ldb, info);
}
inline void ppsv(char* uplo, integer* n, integer* nrhs, dcomplex* ap, dcomplex* b, integer* ldb, integer* info)
{
  zppsv_(uplo, n, nrhs, ap, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline void ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, float* ap, float* afp, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, double* ap, double* afp, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* ap, scomplex* afp, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* ap, dcomplex* afp, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the Cholesky factorization of a real symmetric matrix ---
inline void pptrf(char* uplo, integer* n, float* ap, integer* info)
{
  spptrf_(uplo, n, ap, info);
}
inline void pptrf(char* uplo, integer* n, double* ap, integer* info)
{
  dpptrf_(uplo, n, ap, info);
}
inline void pptrf(char* uplo, integer* n, scomplex* ap, integer* info)
{
  cpptrf_(uplo, n, ap, info);
}
inline void pptrf(char* uplo, integer* n, dcomplex* ap, integer* info)
{
  zpptrf_(uplo, n, ap, info);
}

// --- computes the inverse of a real symmetric matrix ---
inline void pptri(char* uplo, integer* n, float* ap, integer* info)
{
  spptri_(uplo, n, ap, info);
}
inline void pptri(char* uplo, integer* n, double* ap, integer* info)
{
  dpptri_(uplo, n, ap, info);
}
inline void pptri(char* uplo, integer* n, scomplex* ap, integer* info)
{
  cpptri_(uplo, n, ap, info);
}
inline void pptri(char* uplo, integer* n, dcomplex* ap, integer* info)
{
  zpptri_(uplo, n, ap, info);
}

// --- solves a system of linear equations A*X = B with a symmetric matrix ---
inline void pptrs(char* uplo, integer* n, integer* nrhs,  float* ap, float* b, integer* ldb, integer* info)
{
  spptrs_(uplo, n, nrhs, ap, b, ldb, info);
}
inline void pptrs(char* uplo, integer* n, integer* nrhs,  double* ap, double* b, integer* ldb, integer* info)
{
  dpptrs_(uplo, n, nrhs, ap, b, ldb, info);
}
inline void pptrs(char* uplo, integer* n, integer* nrhs,  scomplex* ap, scomplex* b, integer* ldb, integer* info)
{
  cpptrs_(uplo, n, nrhs, ap, b, ldb, info);
}
inline void pptrs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* b, integer* ldb, integer* info)
{
  zpptrs_(uplo, n, nrhs, ap, b, ldb, info);
}

// --- computes the Cholesky factorization ---
inline void pstrf(char* uplo, integer* n, float* a, integer* lda, integer* piv, integer* rank, float* tol, float* work, integer* info)
{
  spstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline void pstrf(char* uplo, integer* n, double* a, integer* lda, integer* piv, integer* rank, double* tol, double* work, integer* info)
{
  dpstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline void pstrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* piv, integer* rank, float* tol, float* work, integer* info)
{
  cpstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline void pstrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* piv, integer* rank, double* tol, double* work, integer* info)
{
  zpstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}

// --- computes the reciprocal of the condition number ---
inline void ptcon(integer* n,  float* d,  float* e, float* anorm, float* rcond, float* rwork, integer* info)
{
  sptcon_(n, d, e, anorm, rcond, rwork, info);
}
inline void ptcon(integer* n,  double* d,  double* e, double* anorm, double* rcond, double* rwork, integer* info)
{
  dptcon_(n, d, e, anorm, rcond, rwork, info);
}
inline void ptcon(integer* n,  float* d,  scomplex* e, float* anorm, float* rcond, float* rwork, integer* info)
{
  cptcon_(n, d, e, anorm, rcond, rwork, info);
}
inline void ptcon(integer* n,  double* d,  dcomplex* e, double* anorm, double* rcond, double* rwork, integer* info)
{
  zptcon_(n, d, e, anorm, rcond, rwork, info);
}

// --- computes all eigenvalues and, optionally, eigenvectors of a symmetric matrix ---
inline void pteqr(char* compz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* info)
{
  spteqr_(compz, n, d, e, z, ldz, work, info);
}
inline void pteqr(char* compz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* info)
{
  dpteqr_(compz, n, d, e, z, ldz, work, info);
}
inline void pteqr(char* compz, integer* n, float* d, float* e, scomplex* z, integer* ldz, float* work, integer* info)
{
  cpteqr_(compz, n, d, e, z, ldz, work, info);
}
inline void pteqr(char* compz, integer* n, double* d, double* e, dcomplex* z, integer* ldz, double* work, integer* info)
{
  zpteqr_(compz, n, d, e, z, ldz, work, info);
}

// --- improves the computed solution to a system of linear equations---
inline void ptrfs(integer* n, integer* nrhs,  float* d,  float* e,  float* df,  float* ef,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* info)
{
  sptrfs_(n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, info);
}
inline void ptrfs(integer* n, integer* nrhs,  double* d,  double* e,  double* df,  double* ef,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* info)
{
  dptrfs_(n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, info);
}
inline void ptrfs(char *uplo, integer* n, integer* nrhs,  float* d,  scomplex* e,  float* df,  scomplex* ef,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cptrfs_(uplo, n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void ptrfs(char *uplo, integer* n, integer* nrhs,  double* d,  dcomplex* e,  double* df,  dcomplex* ef,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zptrfs_(uplo, n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for PT matrices ---
inline void ptsv(integer* n, integer* nrhs, float* d, float* e, float* b, integer* ldb, integer* info)
{
  sptsv_(n, nrhs, d, e, b, ldb, info);
}
inline void ptsv(integer* n, integer* nrhs, double* d, double* e, double* b, integer* ldb, integer* info)
{
  dptsv_(n, nrhs, d, e, b, ldb, info);
}
inline void ptsv(integer* n, integer* nrhs, float* d, scomplex* e, scomplex* b, integer* ldb, integer* info)
{
  cptsv_(n, nrhs, d, e, b, ldb, info);
}
inline void ptsv(integer* n, integer* nrhs, double* d, dcomplex* e, dcomplex* b, integer* ldb, integer* info)
{
  zptsv_(n, nrhs, d, e, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for PT matrices ---
inline void ptsvx(char* fact, integer* n, integer* nrhs,  float* d,  float* e, float* df, float* ef,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* info)
{
  sptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, info);
}
inline void ptsvx(char* fact, integer* n, integer* nrhs,  double* d,  double* e, double* df, double* ef,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* info)
{
  dptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, info);
}
inline void ptsvx(char* fact, integer* n, integer* nrhs,  float* d,  scomplex* e, float* df, scomplex* ef,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void ptsvx(char* fact, integer* n, integer* nrhs,  double* d,  dcomplex* e, double* df, dcomplex* ef,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the L*D*L**T factorization of a real symmetric matrix ---
inline void pttrf(integer* n, float* d, float* e, integer* info)
{
  spttrf_(n, d, e, info);
}
inline void pttrf(integer* n, double* d, double* e, integer* info)
{
  dpttrf_(n, d, e, info);
}
inline void pttrf(integer* n, float* d, scomplex* e, integer* info)
{
  cpttrf_(n, d, e, info);
}
inline void pttrf(integer* n, double* d, dcomplex* e, integer* info)
{
  zpttrf_(n, d, e, info);
}

// --- solves a tridiagonal system of the form A * X = B---
inline void pttrs(integer* n, integer* nrhs,  float* d,  float* e, float* b, integer* ldb, integer* info)
{
  spttrs_(n, nrhs,  d, e, b, ldb, info);
}
inline void pttrs(integer* n, integer* nrhs,  double* d,  double* e, double* b, integer* ldb, integer* info)
{
  dpttrs_(n, nrhs,  d, e, b, ldb, info);
}
inline void pttrs(char *uplo, integer* n, integer* nrhs,  float* d,  scomplex* e, scomplex* b, integer* ldb, integer* info)
{
  cpttrs_(uplo, n, nrhs,  d, e, b, ldb, info);
}
inline void pttrs(char *uplo, integer* n, integer* nrhs,  double* d,  dcomplex* e, dcomplex* b, integer* ldb, integer* info)
{
  zpttrs_(uplo, n, nrhs,  d, e, b, ldb, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void sbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* info)
{
  ssbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info);
}
inline void sbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* info)
{
  dsbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info);
}
inline void hbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  chbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, info);
}
inline void hbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zhbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void sbev(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* info)
{
  ssbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info);
}
inline void sbev(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* info)
{
  dsbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info);
}
inline void hbev(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* info)
{
  chbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info);
}
inline void hbev(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* info)
{
  zhbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void sbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void sbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void hbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  chbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void hbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zhbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void sbevd(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void sbevd(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void hbevd(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  chbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void hbevd(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zhbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void sbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  ssbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline void sbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  dsbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline void hbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, scomplex* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  chbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline void hbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, dcomplex* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zhbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void sbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  ssbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void sbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  dsbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void hbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, scomplex* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  chbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline void hbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, dcomplex* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zhbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- reduces a real symmetric-definite banded generalized eigenproblem  to standard form ---
inline void sbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab,  float* bb, integer* ldbb, float* x, integer* ldx, float* work, integer* info)
{
  ssbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, info);
}
inline void sbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab,  double* bb, integer* ldbb, double* x, integer* ldx, double* work, integer* info)
{
  dsbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, info);
}
inline void hbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab,  scomplex* bb, integer* ldbb, scomplex* x, integer* ldx, scomplex* work, float* rwork, integer* info)
{
  chbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, rwork, info);
}
inline void hbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab,  dcomplex* bb, integer* ldbb, dcomplex* x, integer* ldx, dcomplex* work, double* rwork, integer* info)
{
  zhbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors ---
inline void sbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab, float* bb, integer* ldbb, float* w, float* z, integer* ldz, float* work, integer* info)
{
  ssbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, info);
}
inline void sbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab, double* bb, integer* ldbb, double* w, double* z, integer* ldz, double* work, integer* info)
{
  dsbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, info);
}
inline void hbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab, scomplex* bb, integer* ldbb, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* info)
{
  chbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info);
}
inline void hbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb, integer* ldbb, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* info)
{
  zhbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors ---
inline void sbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab, float* bb, integer* ldbb, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void sbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab, double* bb, integer* ldbb, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void hbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab, scomplex* bb, integer* ldbb, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  chbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void hbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb, integer* ldbb, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zhbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes selected eigenvalues, and optionally, eigenvectors ---
inline void sbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab, float* bb, integer* ldbb, float* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  ssbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void sbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab, double* bb, integer* ldbb, double* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  dsbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void hbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab, scomplex* bb, integer* ldbb, scomplex* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  chbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline void hbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb, integer* ldbb, dcomplex* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zhbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- reduces a real symmetric band matrix A to symmetric tridiagonal form T ---
inline void sbtrd(char* vect, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* d, float* e, float* q, integer* ldq, float* work, integer* info)
{
  ssbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}
inline void sbtrd(char* vect, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* d, double* e, double* q, integer* ldq, double* work, integer* info)
{
  dsbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}
inline void hbtrd(char* vect, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* d, float* e, scomplex* q, integer* ldq, scomplex* work, integer* info)
{
  chbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}
inline void hbtrd(char* vect, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* d, double* e, dcomplex* q, integer* ldq, dcomplex* work, integer* info)
{
  zhbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}

// --- performs a symmetric rank-k operation for matrix in RFP format ---
inline void sfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, float* alpha, float* a, integer* lda, float* beta, float* c)
{
  ssfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}
inline void sfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, double* alpha,  double* a, integer* lda, double* beta, double* c)
{
  dsfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}
inline void hfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, float* alpha,  scomplex* a, integer* lda, float* beta, scomplex* c)
{
  chfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}
inline void hfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, double* alpha,  dcomplex* a, integer* lda, double* beta, dcomplex* c)
{
  zhfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}

// --- estimates the reciprocal of the condition number ---
inline void spcon(char* uplo, integer* n, float* ap, integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  sspcon_(uplo, n, ap, ipiv, anorm, rcond, work, iwork, info);
}
inline void spcon(char* uplo, integer* n, double* ap, integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dspcon_(uplo, n, ap, ipiv, anorm, rcond, work, iwork, info);
}
inline void spcon(char* uplo, integer* n, scomplex* ap, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  cspcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}
inline void spcon(char* uplo, integer* n,  dcomplex* ap, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  zspcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}
inline void hpcon(char* uplo, integer* n, scomplex* ap, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  chpcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}
inline void hpcon(char* uplo, integer* n, dcomplex* ap, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  zhpcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void spev(char *jobz, char *uplo, integer* n, float* ap, float* w, float* z, integer*ldq, float* work, integer* info)
{
  sspev_(jobz, uplo, n, ap, w, z, ldq, work, info);
}
inline void spev(char* jobz, char* uplo, integer* n, double* ap, double* w, double* z, integer* ldq, double* work, integer* info)
{
  dspev_(jobz, uplo, n, ap, w, z, ldq, work, info);
}
inline void hpev(char* jobz, char* uplo, integer* n, scomplex* ap, float* w, scomplex* z, integer* ldq, scomplex* work, float* rwork, integer* info)
{
  chpev_(jobz, uplo, n, ap, w, z, ldq, work, rwork, info);
}
inline void hpev(char* jobz, char* uplo, integer* n, dcomplex* ap, double* w, dcomplex* z, integer* ldq, dcomplex* work, double* rwork, integer* info)
{
  zhpev_(jobz, uplo, n, ap, w, z, ldq, work, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void spevd(char* jobz, char* uplo, integer* n, float* ap, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  sspevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void spevd(char* jobz, char* uplo, integer* n, double* ap, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dspevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void hpevd(char* jobz, char* uplo, integer* n, scomplex* ap, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  chpevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void hpevd(char* jobz, char* uplo, integer* n, dcomplex* ap, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zhpevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void spevx(char* jobz, char* range, char* uplo, integer* n, float* ap, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  sspevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void spevx(char* jobz, char* range, char* uplo, integer* n, double* ap, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  dspevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void hpevx(char* jobz, char* range, char* uplo, integer* n, scomplex* ap, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  chpevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline void hpevx(char* jobz, char* range, char* uplo, integer* n, dcomplex* ap, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zhpevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- reduces a real symmetric-definite generalized eigenproblem to standard form ---
inline void spgst(integer* itype, char* uplo, integer* n, float* ap,  float* bp, integer* info)
{
  sspgst_(itype, uplo, n, ap, bp, info);
}
inline void spgst(integer* itype, char* uplo, integer* n, double* ap,  double* bp, integer* info)
{
  dspgst_(itype, uplo, n, ap, bp, info);
}
inline void hpgst(integer* itype, char* uplo, integer* n, scomplex* ap,  scomplex* bp, integer* info)
{
  chpgst_(itype, uplo, n, ap, bp, info);
}
inline void hpgst(integer* itype, char* uplo, integer* n, dcomplex* ap,  dcomplex* bp, integer* info)
{
  zhpgst_(itype, uplo, n, ap, bp, info);
}

// --- computes all the eigenvalues ---
inline void spgv(integer* itype, char* jobz, char* uplo, integer* n, float* ap, float* bp, float* w, float* z, integer* ldz, float* work, integer* info)
{
  sspgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info);
}
inline void spgv(integer* itype, char* jobz, char* uplo, integer* n, double* ap, double* bp, double* w, double* z, integer* ldz, double* work, integer* info)
{
  dspgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info);
}
inline void hpgv(integer* itype, char* jobz, char* uplo, integer* n, scomplex* ap, scomplex* bp, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* info)
{
  chpgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, info);
}
inline void hpgv(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* ap, dcomplex* bp, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* info)
{
  zhpgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, info);
}

// --- computes all the eigenvalues ---
inline void spgvd(integer* itype, char* jobz, char* uplo, integer* n, float* ap, float* bp, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  sspgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void spgvd(integer* itype, char* jobz, char* uplo, integer* n, double* ap, double* bp, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dspgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline void hpgvd(integer* itype, char* jobz, char* uplo, integer* n, scomplex* ap, scomplex* bp, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  chpgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void hpgvd(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* ap, dcomplex* bp, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zhpgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes all the eigenvalues ---
inline void spgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, float* ap, float* bp, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  sspgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void spgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, double* ap, double* bp, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  dspgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void hpgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, scomplex* ap, scomplex* bp, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  chpgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline void hpgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, dcomplex* ap, dcomplex* bp, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zhpgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline void spsv(char* uplo, integer* n, integer* nrhs, float* ap, integer* ipiv, float* b, integer* ldb, integer* info)
{
  sspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline void spsv(char* uplo, integer* n, integer* nrhs, double* ap, integer* ipiv, double* b, integer* ldb, integer* info)
{
  dspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline void spsv(char* uplo, integer* n, integer* nrhs, scomplex* ap, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  cspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline void spsv(char* uplo, integer* n, integer* nrhs, dcomplex* ap, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline void hpsv(char* uplo, integer* n, integer* nrhs, scomplex* ap, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  chpsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline void hpsv(char* uplo, integer* n, integer* nrhs, dcomplex* ap, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zhpsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline void spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  float* ap, float* afp, integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  sspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  double* ap, double* afp, integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline void spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* ap, scomplex* afp, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* afp, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void hpsvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* ap, scomplex* afp, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  chpsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline void hpsvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* afp, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zhpsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- reduces a real symmetric matrix A stored in packed form to symmetric tridiagonal form T  ---
inline void sptrd(char* uplo, integer* n, float* ap, float* d, float* e, float* tau, integer* info)
{
  ssptrd_(uplo, n, ap, d, e, tau, info);
}
inline void sptrd(char* uplo, integer* n, double* ap, double* d, double* e, double* tau, integer* info)
{
  dsptrd_(uplo, n, ap, d, e, tau, info);
}
inline void hptrd(char* uplo, integer* n, scomplex* ap, float* d, float* e, scomplex* tau, integer* info)
{
  chptrd_(uplo, n, ap, d, e, tau, info);
}
inline void hptrd(char* uplo, integer* n, dcomplex* ap, double* d, double* e, dcomplex* tau, integer* info)
{
  zhptrd_(uplo, n, ap, d, e, tau, info);
}

// --- computes the factorization of a real symmetric matrix  ---
inline void sptrf(char* uplo, integer* n, float* ap, integer* ipiv, integer* info)
{
  ssptrf_(uplo, n, ap, ipiv, info);
}
inline void sptrf(char* uplo, integer* n, double* ap, integer* ipiv, integer* info)
{
  dsptrf_(uplo, n, ap, ipiv, info);
}
inline void sptrf(char* uplo, integer* n, scomplex* ap, integer* ipiv, integer* info)
{
  csptrf_(uplo, n, ap, ipiv, info);
}
inline void sptrf(char* uplo, integer* n, dcomplex* ap, integer* ipiv, integer* info)
{
  zsptrf_(uplo, n, ap, ipiv, info);
}
inline void hptrf(char* uplo, integer* n, scomplex* ap, integer* ipiv, integer* info)
{
  chptrf_(uplo, n, ap, ipiv, info);
}
inline void hptrf(char* uplo, integer* n, dcomplex* ap, integer* ipiv, integer* info)
{
  zhptrf_(uplo, n, ap, ipiv, info);
}

// --- computes the inverse of a real symmetric indefinite matrix ---
inline void sptri(char* uplo, integer* n, float* ap,  integer* ipiv, float* work, integer* info)
{
  ssptri_(uplo, n, ap, ipiv, work, info);
}
inline void sptri(char* uplo, integer* n, double* ap,  integer* ipiv, double* work, integer* info)
{
  dsptri_(uplo, n, ap, ipiv, work, info);
}
inline void sptri(char* uplo, integer* n, scomplex* ap,  integer* ipiv, scomplex* work, integer* info)
{
  csptri_(uplo, n, ap, ipiv, work, info);
}
inline void sptri(char* uplo, integer* n, dcomplex* ap,  integer* ipiv, dcomplex* work, integer* info)
{
  zsptri_(uplo, n, ap, ipiv, work, info);
}
inline void hptri(char* uplo, integer* n, scomplex* ap, integer* ipiv, scomplex* work, integer* info)
{
  chptri_(uplo, n, ap, ipiv, work, info);
}
inline void hptri(char* uplo, integer* n, dcomplex* ap, integer* ipiv, dcomplex* work, integer* info)
{
  zhptri_(uplo, n, ap, ipiv, work, info);
}

// --- solves a system of linear equations A*X = B ---
inline void sptrs(char* uplo, integer* n, integer* nrhs,  float* ap,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  ssptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline void sptrs(char* uplo, integer* n, integer* nrhs,  double* ap,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  dsptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline void sptrs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  csptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline void sptrs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zsptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline void hptrs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  chptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline void hptrs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zhptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}

// --- computes the eigenvalues of a symmetric tridiagonal matrix T ---
inline void stebz(char* range, char* order, integer* n, float* vl, float* vu, integer* il, integer* iu, float* abstol, float* d, float* e, integer* m, integer* nsplit, float* w, integer* iblock, integer* isplit, float* work, integer* iwork, integer* info)
{
  sstebz_(range, order, n, vl, vu, il, iu, abstol, d,  e, m,  nsplit, w, iblock, isplit, work, iwork, info);
}
inline void stebz(char* range, char* order, integer* n, double* vl, double* vu, integer* il, integer* iu, double* abstol, double* d, double* e, integer* m, integer* nsplit, double* w, integer* iblock, integer* isplit, double* work, integer* iwork, integer* info)
{
  dstebz_(range, order, n, vl, vu, il, iu, abstol, d,  e, m,  nsplit, w, iblock, isplit, work, iwork, info);
}

// --- computes selected eigenvalues and, optionally, eigenvectors of a real symmetric tridiagonal matrix T---
inline void stegr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  sstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline void stegr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline void stegr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  cstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline void stegr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  zstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvectors of a real symmetric tridiagonal matrix T ---
inline void stein(integer* n,  float* d,  float* e, integer* m,  float* w,  integer* iblock,  integer* isplit, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  sstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}
inline void stein(integer* n,  double* d,  double* e, integer* m,  double* w,  integer* iblock,  integer* isplit, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  dstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}
inline void stein(integer* n,  float* d,  float* e, integer* m,  float* w,  integer* iblock,  integer* isplit, scomplex* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  cstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}
inline void stein(integer* n,  double* d,  double* e, integer* m,  double* w,  integer* iblock,  integer* isplit, dcomplex* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  zstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}

// --- STERF computes all eigenvalues of a symmetric tridiagonal matrix ---
inline void sterf(integer* n, float* d, float* e, integer* info)
{
  ssterf_(n, d, e, info);
}
inline void sterf(integer* n, double* d, double* e, integer* info)
{
  dsterf_(n, d, e, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void stev(char* jobz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* info)
{
  sstev_(jobz, n, d, e, z, ldz, work, info);
}
inline void stev(char* jobz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* info)
{
  dstev_(jobz, n, d, e, z, ldz, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void stevd(char* jobz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  sstevd_(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}
inline void stevd(char* jobz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dstevd_(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void stevr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  sstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline void stevr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline void stevx(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  sstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline void stevx(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  dstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}

// --- estimates the reciprocal of the condition number of a real symmetric matrix A ---
inline void sycon_3(char* uplo, integer* n,  float* a, integer* lda,  float* e,  integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  ssycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, iwork, info);
}
inline void sycon_3(char* uplo, integer* n,  double* a, integer* lda,  double* e,  integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dsycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, iwork, info);
}
inline void sycon_3(char* uplo, integer* n,  scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  csycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, info);
}
inline void sycon_3(char* uplo, integer* n,  dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  zsycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, info);
}

// --- estimates the reciprocal of the condition number of a real symmetric matrix A ---
inline void sycon(char* uplo, integer* n,  float* a, integer* lda,  integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  ssycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, iwork, info);
}
inline void sycon(char* uplo, integer* n,  double* a, integer* lda,  integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dsycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, iwork, info);
}
inline void sycon(char* uplo, integer* n,  scomplex* a, integer* lda,  integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  csycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, info);
}
inline void sycon(char* uplo, integer* n,  dcomplex* a, integer* lda,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  zsycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, info);
}

// --- convert A given by TRF into L and D and vice-versa ---
inline void syconv(char* uplo, char* way, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* info)
{
  ssyconv_(uplo, way, n, a, lda, ipiv, work, info);
}
inline void syconv(char* uplo, char* way, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* info)
{
  dsyconv_(uplo, way, n, a, lda, ipiv, work, info);
}
inline void syconv(char* uplo, char* way, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* info)
{
  csyconv_(uplo, way, n, a, lda, ipiv, work, info);
}
inline void syconv(char* uplo, char* way, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* info)
{
  zsyconv_(uplo, way, n, a, lda, ipiv, work, info);
}

// --- computes row and column scalings intended to equilibrate a symmetric matrix A ---
inline void syequb(char* uplo, integer* n,  float* a, integer* lda, float* s, float* scond, float* amax, float* work, integer* info)
{
  ssyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline void syequb(char* uplo, integer* n,  double* a, integer* lda, double* s, double* scond, double* amax, double* work, integer* info)
{
  dsyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline void syequb(char* uplo, integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, scomplex* work, integer* info)
{
  csyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline void syequb(char* uplo, integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, dcomplex* work, integer* info)
{
  zsyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline void syev_2stage(char* jobz, char* uplo, integer* n, float* a, integer* lda, float* w, float* work, integer* lwork, integer* info)
{
  ssyev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, info);
}
inline void syev_2stage(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, integer* info)
{
  dsyev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline void syevd_2stage(char* jobz, char* uplo, integer* n, float* a, integer* lda, float* w, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssyevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}
inline void syevd_2stage(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsyevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices --
inline void syevr_2stage(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssyevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline void syevr_2stage(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsyevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline void syevx_2stage(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  ssyevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline void syevx_2stage(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  dsyevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline void syevx(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  ssyevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline void syevx(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  dsyevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}

// --- computes all the eigenvalues, the eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline void sygv_2stage(integer* itype, char* jobz, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* w, float* work, integer* lwork, integer* info)
{
  ssygv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}
inline void sygv_2stage(integer* itype, char* jobz, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* w, double* work, integer* lwork, integer* info)
{
  dsygv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}

// --- computes all the eigenvalues, the eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline void sygv(integer* itype, char* jobz, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* w, float* work, integer* lwork, integer* info)
{
  ssygv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}
inline void sygv(integer* itype, char* jobz, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* w, double* work, integer* lwork, integer* info)
{
  dsygv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}

// --- computes all the eigenvalues, the eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline void sygvd(integer* itype, char* jobz, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* w, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ssygvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info);
}
inline void sygvd(integer* itype, char* jobz, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* w, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dsygvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info);
}

// --- computes selected eigenvalues, eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline void sygvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  ssygvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline void sygvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  dsygvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}

// --- improves the computed solution to a system of linear equations---
inline void syrfs(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  ssyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void syrfs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dsyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void syrfs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  csyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void syrfs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zsyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- improves the computed solution to a system of linear equation---
inline void syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* s,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  ssyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* s,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dsyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  float* s,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  csyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  double* s,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zsyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline void sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* tb, integer* ltb, integer* ipiv, integer* ipiv2, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  ssysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline void sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* tb, integer* ltb, integer* ipiv, integer* ipiv2, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dsysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline void sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  csysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline void sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zsysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline void sysv_aa(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  ssysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_aa(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dsysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_aa(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  csysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_aa(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zsysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline void sysv_rk(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* e, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  ssysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_rk(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* e, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dsysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_rk(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  csysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_rk(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zsysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline void sysv_rook(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  ssysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_rook(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dsysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_rook(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  csysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv_rook(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zsysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline void sysv(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  ssysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dsysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  csysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void sysv(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zsysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline void sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  float* a, integer* lda, float* af, integer* ldaf, integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* lwork, integer* iwork, integer* info)
{
  ssysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, iwork, info);
}
inline void sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  double* a, integer* lda, double* af, integer* ldaf, integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* lwork, integer* iwork, integer* info)
{
  dsysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, iwork, info);
}
inline void sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  csysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}
inline void sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zsysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}

// --- uses the diagonal pivoting factorization to compute the solution to a real system of linear equations ---
inline void sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, integer* ipiv, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  ssysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, integer* ipiv, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  dsysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline void sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  csysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zsysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- applies an elementary permutation on the rows and columns of a symmetric matrix ---
inline void syswapr(char* uplo, integer* n, float* a, integer* lda, integer* i1, integer* i2)
{
  ssyswapr_(uplo, n, a, lda, i1, i2);
}
inline void syswapr(char* uplo, integer* n, double* a, integer* lda, integer* i1, integer* i2)
{
  dsyswapr_(uplo, n, a, lda, i1, i2);
}
inline void syswapr(char* uplo, integer* n, scomplex* a, integer* lda, integer* i1, integer* i2)
{
  csyswapr_(uplo, n, a, lda, i1, i2);
}
inline void syswapr(char* uplo, integer* n, dcomplex* a, integer* lda, integer* i1, integer* i2)
{
  zsyswapr_(uplo, n, a, lda, i1, i2);
}

// --- computes the factorization of a real symmetric matrix A using the Aasen's algorithm---
inline void sytrf_aa_2stage(char* uplo, integer* n, float* a, integer* lda, float* tb, integer* ltb, integer* ipiv, integer* ipiv2, float* work, integer* lwork, integer* info)
{
  ssytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline void sytrf_aa_2stage(char* uplo, integer* n, double* a, integer* lda, double* tb, integer* ltb, integer* ipiv, integer* ipiv2, double* work, integer* lwork, integer* info)
{
  dsytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline void sytrf_aa_2stage(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* work, integer* lwork, integer* info)
{
  csytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline void sytrf_aa_2stage(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* work, integer* lwork, integer* info)
{
  zsytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the Aasen's algorithm ---
inline void sytrf_aa(char* uplo, integer* n, float* a, integer* lda, integer* ipiv, float* work, integer* lwork, integer* info)
{
  ssytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf_aa(char* uplo, integer* n, double* a, integer* lda, integer* ipiv, double* work, integer* lwork, integer* info)
{
  dsytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf_aa(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  csytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf_aa(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zsytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline void sytrf_rk(char* uplo, integer* n, float* a, integer* lda, float* e, integer* ipiv, float* work, integer* lwork, integer* info)
{
  ssytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline void sytrf_rk(char* uplo, integer* n, double* a, integer* lda, double* e, integer* ipiv, double* work, integer* lwork, integer* info)
{
  dsytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline void sytrf_rk(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  csytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline void sytrf_rk(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zsytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting method---
inline void sytrf_rook(char* uplo, integer* n, float* a, integer* lda, integer* ipiv, float* work, integer* lwork, integer* info)
{
  ssytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf_rook(char* uplo, integer* n, double* a, integer* lda, integer* ipiv, double* work, integer* lwork, integer* info)
{
  dsytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf_rook(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  csytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf_rook(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zsytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method ---
inline void sytrf(char* uplo, integer* n, float* a, integer* lda, integer* ipiv, float* work, integer* lwork, integer* info)
{
  ssytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf(char* uplo, integer* n, double* a, integer* lda, integer* ipiv, double* work, integer* lwork, integer* info)
{
  dsytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  csytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zsytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- SYTRI_3 computes the inverse of a real symmetric indefinite matrix A using the factorization computed by SYTRF_RK or SYTRF_BK---
inline void sytri_3(char* uplo, integer* n, float* a, integer* lda,  float* e,  integer* ipiv, float* work, integer* lwork, integer* info)
{
  ssytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}
inline void sytri_3(char* uplo, integer* n, double* a, integer* lda,  double* e,  integer* ipiv, double* work, integer* lwork, integer* info)
{
  dsytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}
inline void sytri_3(char* uplo, integer* n, scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  csytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}
inline void sytri_3(char* uplo, integer* n, dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zsytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}

// --- computes the inverse of a real symmetric indefinite matrix A---
inline void sytri(char* uplo, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* info)
{
  ssytri_(uplo, n, a, lda, ipiv, work, info);
}
inline void sytri(char* uplo, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* info)
{
  dsytri_(uplo, n, a, lda, ipiv, work, info);
}
inline void sytri(char* uplo, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* info)
{
  csytri_(uplo, n, a, lda, ipiv, work, info);
}
inline void sytri(char* uplo, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* info)
{
  zsytri_(uplo, n, a, lda, ipiv, work, info);
}

// --- computes the inverse of a REAL symmetric indefinite matrix A using the factorization ---
inline void sytri2(char* uplo, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* lwork, integer* info)
{
  ssytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytri2(char* uplo, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* lwork, integer* info)
{
  dsytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytri2(char* uplo, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  csytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void sytri2(char* uplo, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zsytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the inverse of a real symmetric indefinite matrix A using the factorization ---
inline void sytri2x(char* uplo, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* nb, integer* info)
{
  ssytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline void sytri2x(char* uplo, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* nb, integer* info)
{
  dsytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline void sytri2x(char* uplo, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* nb, integer* info)
{
  csytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline void sytri2x(char* uplo, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* nb, integer* info)
{
  zsytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}

// --- solves a system of linear equations A * X = B with a real symmetric matrix A ---
inline void sytrs_3(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  float* e,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  ssytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}
inline void sytrs_3(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  double* e,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  dsytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}
inline void sytrs_3(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  csytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}
inline void sytrs_3(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zsytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A---
inline void sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* tb, integer* ltb, integer* ipiv, integer* ipiv2, float* b, integer* ldb, integer* info)
{
  ssytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline void sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* tb, integer* ltb, integer* ipiv, integer* ipiv2, double* b, integer* ldb, integer* info)
{
  dsytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline void sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, integer* info)
{
  csytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline void sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, integer* info)
{
  zsytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A ---
inline void sytrs_aa(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  ssytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline void sytrs_aa(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  dsytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline void sytrs_aa(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  csytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline void sytrs_aa(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zsytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A---
inline void sytrs_rook(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  ssytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void sytrs_rook(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  dsytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void sytrs_rook(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  csytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void sytrs_rook(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zsytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A ---
inline void sytrs(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  ssytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void sytrs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  dsytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void sytrs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  csytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void sytrs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zsytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A ---
inline void sytrs2(char* uplo, integer* n, integer* nrhs, float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, float* work, integer* info)
{
  ssytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}
inline void sytrs2(char* uplo, integer* n, integer* nrhs, double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, double* work, integer* info)
{
  dsytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}
inline void sytrs2(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* info)
{
  csytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}
inline void sytrs2(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* info)
{
  zsytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}

// --- computes the eigenvalues of a real matrix pair (H,T)---
inline void hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, float* h, integer* ldh, float* t, integer* ldt, float* alphar, float* alphai, float* beta, float* q, integer* ldq, float* z, integer* ldz, float* work, integer* lwork, integer* info)
{
  shgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
}
inline void hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, double* h, integer* ldh, double* t, integer* ldt, double* alphar, double* alphai, double* beta, double* q, integer* ldq, double* z, integer* ldz, double* work, integer* lwork, integer* info)
{
  dhgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
}
inline void hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* h, integer* ldh, scomplex* t, integer* ldt, scomplex* alpha, scomplex* beta, scomplex* q, integer* ldq, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  chgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
}
inline void hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* h, integer* ldh, dcomplex* t, integer* ldt, dcomplex* alpha, dcomplex* beta, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zhgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
}

// --- uses inverse iteration to find specified right and/or left eigenvectors of a real upper Hessenberg matrix H ---
inline void hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  float* h, integer* ldh, float* wr,  float* wi, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* mm, integer* m, float* work, integer* ifaill, integer* ifailr, integer* info)
{
  shsein_(job, eigsrc, initv, select, n,  h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work, ifaill, ifailr, info);
}
inline void hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  double* h, integer* ldh, double* wr,  double* wi, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* mm, integer* m, double* work, integer* ifaill, integer* ifailr, integer* info)
{
  dhsein_(job, eigsrc, initv, select, n,  h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work, ifaill, ifailr, info);
}
inline void hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  scomplex* h, integer* ldh, scomplex* w, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* mm, integer* m, scomplex* work, float* rwork, integer* ifaill, integer* ifailr, integer* info)
{
  chsein_(job, eigsrc, initv, select, n,  h, ldh, w,  vl, ldvl, vr, ldvr, mm, m, work, rwork, ifaill, ifailr, info);
}
inline void hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  dcomplex* h, integer* ldh, dcomplex* w, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* mm, integer* m, dcomplex* work, double* rwork, integer* ifaill, integer* ifailr, integer* info)
{
  zhsein_(job, eigsrc, initv, select, n,  h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, rwork, ifaill, ifailr, info);
}

// --- computes the eigenvalues of a Hessenberg matrix H ---
inline void hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, float* h, integer* ldh, float* wr, float* wi, float *z, integer* ldz, float* work, integer* lwork, integer* info)
{
  shseqr_(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
}
inline void hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, double* h, integer* ldh, double* wr, double* wi, double *z, integer* ldz, double* work, integer* lwork, integer* info)
{
  dhseqr_(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
}
inline void hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* h, integer* ldh, scomplex* w, scomplex*z, integer* ldz, scomplex* work, integer* lwork, integer* info)
{
  chseqr_(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
}
inline void hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* h, integer* ldh, dcomplex* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, integer* info)
{
  zhseqr_(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
}

// --- TBCON estimates the reciprocal of the condition number of a triangular band matrix A---
inline void tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  float* ab, integer* ldab, float* rcond, float* work, integer* iwork, integer* info)
{
  stbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, iwork, info);
}
inline void tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  double* ab, integer* ldab, double* rcond, double* work, integer* iwork, integer* info)
{
  dtbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, iwork, info);
}
inline void tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  scomplex* ab, integer* ldab, float* rcond, scomplex* work, float* rwork, integer* info)
{
  ctbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, rwork, info);
}
inline void tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  dcomplex* ab, integer* ldab, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  ztbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, rwork, info);
}

// --- TBRFS provides error bounds and backward error estimates for the solution to a system of linear equations ---
inline void tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab,  float* b, integer* ldb,  float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  stbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline void tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab,  double* b, integer* ldb,  double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dtbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline void tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab,  scomplex* b, integer* ldb,  scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  ctbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}
inline void tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab,  dcomplex* b, integer* ldb,  dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  ztbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}

// --- solves a triangular system of the form A * X = B  or  A**T * X = B ---
inline void tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab, float* b, integer* ldb, integer* info)
{
  stbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}
inline void tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab, double* b, integer* ldb, integer* info)
{
  dtbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}
inline void tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab, scomplex* b, integer* ldb, integer* info)
{
  ctbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}
inline void tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab, dcomplex* b, integer* ldb, integer* info)
{
  ztbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}

// --- solves a matrix equation (one operand is a triangular matrix in RFP format). ---
inline void tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, float* alpha, float* a, float* b, integer* ldb)
{
  stfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}
inline void tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, double* alpha, double* a, double* b, integer* ldb)
{
  dtfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}
inline void tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, scomplex* alpha, scomplex* a, scomplex* b, integer* ldb)
{
  ctfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}
inline void tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, dcomplex* alpha, dcomplex* a, dcomplex* b, integer* ldb)
{
  ztfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}

// --- computes the inverse of a triangular matrix A stored in RFP format ---
inline void tftri(char* transr, char* uplo, char* diag, integer* n, float* a, integer* info)
{
  stftri_(transr, uplo, diag, n, a, info);
}
inline void tftri(char* transr, char* uplo, char* diag, integer* n, double* a, integer* info)
{
  dtftri_(transr, uplo, diag, n, a, info);
}
inline void tftri(char* transr, char* uplo, char* diag, integer* n, scomplex* a, integer* info)
{
  ctftri_(transr, uplo, diag, n, a, info);
}
inline void tftri(char* transr, char* uplo, char* diag, integer* n, dcomplex* a, integer* info)
{
  ztftri_(transr, uplo, diag, n, a, info);
}

// --- copies a triangular matrix from the rectangular full packed format (TF) to the standard packed format (TP) ---
inline void tfttp(char* transr, char* uplo, integer* n,  float* arf, float* ap, integer* info)
{
  stfttp_(transr, uplo, n,  arf, ap, info);
}
inline void tfttp(char* transr, char* uplo, integer* n,  double* arf, double* ap, integer* info)
{
  dtfttp_(transr, uplo, n,  arf, ap, info);
}
inline void tfttp(char* transr, char* uplo, integer* n,  scomplex* arf, scomplex* ap, integer* info)
{
  ctfttp_(transr, uplo, n,  arf, ap, info);
}
inline void tfttp(char* transr, char* uplo, integer* n,  dcomplex* arf, dcomplex* ap, integer* info)
{
  ztfttp_(transr, uplo, n,  arf, ap, info);
}

// --- copies a triangular matrix from the rectangular full packed format (TF) to the standard full format (TR) ---
inline void tfttr(char* transr, char* uplo, integer* n,  float* arf, float* a, integer* lda, integer* info)
{
  stfttr_(transr, uplo, n,  arf, a, lda, info);
}
inline void tfttr(char* transr, char* uplo, integer* n,  double* arf, double* a, integer* lda, integer* info)
{
  dtfttr_(transr, uplo, n,  arf, a, lda, info);
}
inline void tfttr(char* transr, char* uplo, integer* n,  scomplex* arf, scomplex* a, integer* lda, integer* info)
{
  ctfttr_(transr, uplo, n,  arf, a, lda, info);
}
inline void tfttr(char* transr, char* uplo, integer* n,  dcomplex* arf, dcomplex* a, integer* lda, integer* info)
{
  ztfttr_(transr, uplo, n,  arf, a, lda, info);
}

// --- computes some or all of the right and/or left eigenvectors of a pair of real matrices ---
inline void tgevc(char* side, char* howmny,  logical* select, integer* n,  float* s, integer* lds,  float* p, integer* ldp, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* mm, integer* m, float* work, integer* info)
{
  stgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline void tgevc(char* side, char* howmny,  logical* select, integer* n,  double* s, integer* lds,  double* p, integer* ldp, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* mm, integer* m, double* work, integer* info)
{
  dtgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline void tgevc(char* side, char* howmny,  logical* select, integer* n,  scomplex* s, integer* lds,  scomplex* p, integer* ldp, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* mm, integer* m, scomplex* work, float* rwork, integer* info)
{
  ctgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}
inline void tgevc(char* side, char* howmny,  logical* select, integer* n,  dcomplex* s, integer* lds,  dcomplex* p, integer* ldp, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* mm, integer* m, dcomplex* work, double* rwork, integer* info)
{
  ztgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}

// --- reorders the generalized real Schur decomposition of a real matrix pair ---
inline void tgexc(logical* wantq, logical* wantz, integer* n, float* a, integer* lda, float* b, integer* ldb, float* q, integer* ldq, float* z, integer* ldz, integer* ifst, integer* ilst, float* work, integer* lwork, integer* info)
{
  stgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info);
}
inline void tgexc(logical* wantq, logical* wantz, integer* n, double* a, integer* lda, double* b, integer* ldb, double* q, integer* ldq, double* z, integer* ldz, integer* ifst, integer* ilst, double* work, integer* lwork, integer* info)
{
  dtgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info);
}
inline void tgexc(logical* wantq, logical* wantz, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* q, integer* ldq, scomplex* z, integer* ldz, integer* ifst, integer* ilst, integer* info)
{
  ctgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info);
}
inline void tgexc(logical* wantq, logical* wantz, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, integer* ifst, integer* ilst, integer* info)
{
  ztgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info);
}

// --- reorders the generalized real Schur decomposition of a real matrix pair ---
inline void tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* q, integer* ldq, float* z, integer* ldz, integer* m, float* pl, float* pr, float* dif, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  stgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}
inline void tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* q, integer* ldq, double* z, integer* ldz, integer* m, double* pl, double* pr, double* dif, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dtgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}
inline void tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* q, integer* ldq, scomplex* z, integer* ldz, integer* m, float* pl, float* pr, float* dif, scomplex* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ctgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}
inline void tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, integer* m, double* pl, double* pr, double* dif, dcomplex* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  ztgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}

// --- computes the generalized singular value decomposition (GSVD) of two real upper triangular ---
inline void tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, float* a, integer* lda, float* b, integer* ldb, float* tola, float* tolb, float* alpha, float* beta, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, float* work, integer* ncycle, integer* info)
{
  stgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}
inline void tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, double* a, integer* lda, double* b, integer* ldb, double* tola, double* tolb, double* alpha, double* beta, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, double* work, integer* ncycle, integer* info)
{
  dtgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}
inline void tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* tola, float* tolb, float* alpha, float* beta, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, scomplex* work, integer* ncycle, integer* info)
{
  ctgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}
inline void tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* tola, double* tolb, double* alpha, double* beta, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, dcomplex* work, integer* ncycle, integer* info)
{
  ztgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}

// --- estimates reciprocal condition numbers for specified eigenvalues and/or eigenvectors of a matrix pair ---
inline void tgsna(char* job, char* howmny,  logical* select, integer* n,  float* a, integer* lda,  float* b, integer* ldb,  float* vl, integer* ldvl,  float* vr, integer* ldvr, float* s, float* dif, integer* mm, integer* m, float* work, integer* lwork, integer* iwork, integer* info)
{
  stgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}
inline void tgsna(char* job, char* howmny,  logical* select, integer* n,  double* a, integer* lda,  double* b, integer* ldb,  double* vl, integer* ldvl,  double* vr, integer* ldvr, double* s, double* dif, integer* mm, integer* m, double* work, integer* lwork, integer* iwork, integer* info)
{
  dtgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}
inline void tgsna(char* job, char* howmny,  logical* select, integer* n,  scomplex* a, integer* lda,  scomplex* b, integer* ldb,  scomplex* vl, integer* ldvl,  scomplex* vr, integer* ldvr, float* s, float* dif, integer* mm, integer* m, scomplex* work, integer* lwork, integer* iwork, integer* info)
{
  ctgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}
inline void tgsna(char* job, char* howmny,  logical* select, integer* n,  dcomplex* a, integer* lda,  dcomplex* b, integer* ldb,  dcomplex* vl, integer* ldvl,  dcomplex* vr, integer* ldvr, double* s, double* dif, integer* mm, integer* m, dcomplex* work, integer* lwork, integer* iwork, integer* info)
{
  ztgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}

// ---  solves the generalized Sylvester equation ---
inline void tgsyl(char* trans, integer* ijob, integer* m, integer* n,  float* a, integer* lda,  float* b, integer* ldb, float* c, integer* ldc,  float* d, integer* ldd,  float* e, integer* lde, float* f, integer* ldf, float* scale, float* dif, float* work, integer* lwork, integer* iwork, integer* info)
{
  stgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}
inline void tgsyl(char* trans, integer* ijob, integer* m, integer* n,  double* a, integer* lda,  double* b, integer* ldb, double* c, integer* ldc,  double* d, integer* ldd,  double* e, integer* lde, double* f, integer* ldf, double* scale, double* dif, double* work, integer* lwork, integer* iwork, integer* info)
{
  dtgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}
inline void tgsyl(char* trans, integer* ijob, integer* m, integer* n,  scomplex* a, integer* lda,  scomplex* b, integer* ldb, scomplex* c, integer* ldc,  scomplex* d, integer* ldd,  scomplex* e, integer* lde, scomplex* f, integer* ldf, float* scale, float* dif, scomplex* work, integer* lwork, integer* iwork, integer* info)
{
  ctgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}
inline void tgsyl(char* trans, integer* ijob, integer* m, integer* n,  dcomplex* a, integer* lda,  dcomplex* b, integer* ldb, dcomplex* c, integer* ldc,  dcomplex* d, integer* ldd,  dcomplex* e, integer* lde, dcomplex* f, integer* ldf, double* scale, double* dif, dcomplex* work, integer* lwork, integer* iwork, integer* info)
{
  ztgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}

// --- estimates the reciprocal of the condition number of a packed triangular matrix A ---
inline void tpcon(char* norm, char* uplo, char* diag, integer* n,  float* ap, float* rcond, float* work, integer* iwork, integer* info)
{
  stpcon_(norm, uplo, diag, n,  ap, rcond, work, iwork, info);
}
inline void tpcon(char* norm, char* uplo, char* diag, integer* n,  double* ap, double* rcond, double* work, integer* iwork, integer* info)
{
  dtpcon_(norm, uplo, diag, n,  ap, rcond, work, iwork, info);
}
inline void tpcon(char* norm, char* uplo, char* diag, integer* n,  scomplex* ap, float* rcond, scomplex* work, float* rwork, integer* info)
{
  ctpcon_(norm, uplo, diag, n,  ap, rcond, work, rwork, info);
}
inline void tpcon(char* norm, char* uplo, char* diag, integer* n,  dcomplex* ap, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  ztpcon_(norm, uplo, diag, n,  ap, rcond, work, rwork, info);
}

// --- TPMQRT applies a real orthogonal matrix Q obtained from a triangular-pentagonal ---
inline void tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  float* v, integer* ldv,  float* t, integer* ldt, float* a, integer* lda, float* b, integer* ldb, float* work, integer* info)
{
  stpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline void tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  double* v, integer* ldv,  double* t, integer* ldt, double* a, integer* lda, double* b, integer* ldb, double* work, integer* info)
{
  dtpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline void tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* info)
{
  ctpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline void tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* info)
{
  ztpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}

// --- TPQRT computes a blocked QR factorization of a real triangular-pentagonal matrix ---
inline void tpqrt(integer* m, integer* n, integer* l, integer* nb, float* a, integer* lda, float* b, integer* ldb, float* t, integer* ldt, float* work, integer* info)
{
  stpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}
inline void tpqrt(integer* m, integer* n, integer* l, integer* nb, double* a, integer* lda, double* b, integer* ldb, double* t, integer* ldt, double* work, integer* info)
{
  dtpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}
inline void tpqrt(integer* m, integer* n, integer* l, integer* nb, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* t, integer* ldt, scomplex* work, integer* info)
{
  ctpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}
inline void tpqrt(integer* m, integer* n, integer* l, integer* nb, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* t, integer* ldt, dcomplex* work, integer* info)
{
  ztpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}

// --- computes a QR factorization of a real or complex "triangular-pentagonal" matrix ---
inline void tpqrt2(integer* m, integer* n, integer* l, float* a, integer* lda, float* b, integer* ldb, float* t, integer* ldt, integer* info)
{
  stpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline void tpqrt2(integer* m, integer* n, integer* l, double* a, integer* lda, double* b, integer* ldb, double* t, integer* ldt, integer* info)
{
  dtpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline void tpqrt2(integer* m, integer* n, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* t, integer* ldt, integer* info)
{
  ctpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline void tpqrt2(integer* m, integer* n, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* t, integer* ldt, integer* info)
{
  ztpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}

// --- applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex matrix ---
inline void tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  float* v, integer* ldv,  float* t, integer* ldt, float* a, integer* lda, float* b, integer* ldb, float* work, integer* ldwork)
{
  stprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline void tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  double* v, integer* ldv,  double* t, integer* ldt, double* a, integer* lda, double* b, integer* ldb, double* work, integer* ldwork)
{
  dtprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline void tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* ldwork)
{
  ctprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline void tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* ldwork)
{
  ztprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}

// --- provides error bounds and backward error estimates for the solution to a system of linear equations ---
inline void tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* ap,  float* b, integer* ldb,  float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  stprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline void tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* ap,  double* b, integer* ldb,  double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dtprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline void tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* ap,  scomplex* b, integer* ldb,  scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  ctprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}
inline void tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* b, integer* ldb,  dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  ztprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the inverse of a real upper or lower triangular matrix A stored in packed format ---
inline void tptri(char* uplo, char* diag, integer* n, float* ap, integer* info)
{
  stptri_(uplo, diag, n, ap, info);
}
inline void tptri(char* uplo, char* diag, integer* n, double* ap, integer* info)
{
  dtptri_(uplo, diag, n, ap, info);
}
inline void tptri(char* uplo, char* diag, integer* n, scomplex* ap, integer* info)
{
  ctptri_(uplo, diag, n, ap, info);
}
inline void tptri(char* uplo, char* diag, integer* n, dcomplex* ap, integer* info)
{
  ztptri_(uplo, diag, n, ap, info);
}

// --- solves a triangular system of the form A * X = B  or  A**T * X = B ---
inline void tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* ap, float* b, integer* ldb, integer* info)
{
  stptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}
inline void tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* ap, double* b, integer* ldb, integer* info)
{
  dtptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}
inline void tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* ap, scomplex* b, integer* ldb, integer* info)
{
  ctptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}
inline void tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* b, integer* ldb, integer* info)
{
  ztptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}

// --- copies a triangular matrix from the standard packed format (TP) to the rectangular full packed format (TF) ---
inline void tpttf(char* transr, char* uplo, integer* n,  float* ap, float* arf, integer* info)
{
  stpttf_(transr, uplo, n,  ap, arf, info);
}
inline void tpttf(char* transr, char* uplo, integer* n,  double* ap, double* arf, integer* info)
{
  dtpttf_(transr, uplo, n,  ap, arf, info);
}
inline void tpttf(char* transr, char* uplo, integer* n,  scomplex* ap, scomplex* arf, integer* info)
{
  ctpttf_(transr, uplo, n,  ap, arf, info);
}
inline void tpttf(char* transr, char* uplo, integer* n,  dcomplex* ap, dcomplex* arf, integer* info)
{
  ztpttf_(transr, uplo, n,  ap, arf, info);
}

// --- estimates the reciprocal of the condition number of a triangular matrix A ---
inline void trcon(char* norm, char* uplo, char* diag, integer* n,  float* a, integer* lda, float* rcond, float* work, integer* iwork, integer* info)
{
  strcon_(norm, uplo, diag, n,  a, lda, rcond, work, iwork, info);
}
inline void trcon(char* norm, char* uplo, char* diag, integer* n,  double* a, integer* lda, double* rcond, double* work, integer* iwork, integer* info)
{
  dtrcon_(norm, uplo, diag, n,  a, lda, rcond, work, iwork, info);
}
inline void trcon(char* norm, char* uplo, char* diag, integer* n,  scomplex* a, integer* lda, float* rcond, scomplex* work, float* rwork, integer* info)
{
  ctrcon_(norm, uplo, diag, n,  a, lda, rcond, work, rwork, info);
}
inline void trcon(char* norm, char* uplo, char* diag, integer* n,  dcomplex* a, integer* lda, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  ztrcon_(norm, uplo, diag, n,  a, lda, rcond, work, rwork, info);
}

// --- computes some or all of the right and/or left eigenvectors of a real upper quasi-triangular matrix T ---
inline void trevc(char* side, char* howmny, logical* select, integer* n,  float* t, integer* ldt, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* mm, integer* m, float* work, integer* info)
{
  strevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline void trevc(char* side, char* howmny, logical* select, integer* n,  double* t, integer* ldt, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* mm, integer* m, double* work, integer* info)
{
  dtrevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline void trevc(char* side, char* howmny, logical* select, integer* n,  scomplex* t, integer* ldt, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* mm, integer* m, scomplex* work, float* rwork, integer* info)
{
  ctrevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}
inline void trevc(char* side, char* howmny, logical* select, integer* n,  dcomplex* t, integer* ldt, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* mm, integer* m, dcomplex* work, double* rwork, integer* info)
{
  ztrevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}

// --- reorders the real Schur factorization of a real matrix A = Q*T*Q**T ---
inline void trexc(char* compq, integer* n, float* t, integer* ldt, float* q, integer* ldq, integer* ifst, integer* ilst, float* work, integer* info)
{
  strexc_(compq, n, t, ldt, q, ldq, ifst, ilst, work, info);
}
inline void trexc(char* compq, integer* n, double* t, integer* ldt, double* q, integer* ldq, integer* ifst, integer* ilst, double* work, integer* info)
{
  dtrexc_(compq, n, t, ldt, q, ldq, ifst, ilst, work, info);
}
inline void trexc(char* compq, integer* n, scomplex* t, integer* ldt, scomplex* q, integer* ldq, integer* ifst, integer* ilst, integer* info)
{
  ctrexc_(compq, n, t, ldt, q, ldq, ifst, ilst, info);
}
inline void trexc(char* compq, integer* n, dcomplex* t, integer* ldt, dcomplex* q, integer* ldq, integer* ifst, integer* ilst, integer* info)
{
  ztrexc_(compq, n, t, ldt, q, ldq, ifst, ilst, info);
}

// --- provides error bounds and backward error estimates for the solution to a system of linear equations ---
inline void trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* a, integer* lda,  float* b, integer* ldb,  float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  strrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline void trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* a, integer* lda,  double* b, integer* ldb,  double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dtrrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline void trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* b, integer* ldb,  scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  ctrrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}
inline void trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* b, integer* ldb,  dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  ztrrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}

// --- reorders the real Schur factorization of a real matrix A = Q*T*Q**T ---
inline void trsen(char* job, char* compq,  logical* select, integer* n, float* t, integer* ldt, float* q, integer* ldq, float* wr, float* wi, integer* m, float* s, float* sep, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  strsen_(job, compq,  select, n, t, ldt, q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info);
}
inline void trsen(char* job, char* compq,  logical* select, integer* n, double* t, integer* ldt, double* q, integer* ldq, double* wr, double* wi, integer* m, double* s, double* sep, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  dtrsen_(job, compq,  select, n, t, ldt, q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info);
}
inline void trsen(char* job, char* compq,  logical* select, integer* n, scomplex* t, integer* ldt, scomplex* q, integer* ldq, scomplex* w, integer* m, float* s, float* sep, scomplex* work, integer* lwork, integer* info)
{
  ctrsen_(job, compq,  select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info);
}
inline void trsen(char* job, char* compq,  logical* select, integer* n, dcomplex* t, integer* ldt, dcomplex* q, integer* ldq, dcomplex* w, integer* m, double* s, double* sep, dcomplex* work, integer* lwork, integer* info)
{
  ztrsen_(job, compq,  select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info);
}

// --- estimates reciprocal condition numbers for specified eigenvalues ---
inline void trsna(char* job, char* howmny,  logical* select, integer* n,  float* t, integer* ldt,  float* vl, integer* ldvl,  float* vr, integer* ldvr, float* s, float* sep, integer* mm, integer* m, float* work, integer* ldwork, integer* iwork, integer* info)
{
  strsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, iwork, info);
}
inline void trsna(char* job, char* howmny,  logical* select, integer* n,  double* t, integer* ldt,  double* vl, integer* ldvl,  double* vr, integer* ldvr, double* s, double* sep, integer* mm, integer* m, double* work, integer* ldwork, integer* iwork, integer* info)
{
  dtrsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, iwork, info);
}
inline void trsna(char* job, char* howmny,  logical* select, integer* n,  scomplex* t, integer* ldt,  scomplex* vl, integer* ldvl,  scomplex* vr, integer* ldvr, float* s, float* sep, integer* mm, integer* m, scomplex* work, integer* ldwork, float* rwork, integer* info)
{
  ctrsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, rwork, info);
}
inline void trsna(char* job, char* howmny,  logical* select, integer* n,  dcomplex* t, integer* ldt,  dcomplex* vl, integer* ldvl,  dcomplex* vr, integer* ldvr, double* s, double* sep, integer* mm, integer* m, dcomplex* work, integer* ldwork, double* rwork, integer* info)
{
  ztrsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, rwork, info);
}

// --- solves a triangular system of the form A * X = B  or  A**T * X = B ---
inline void trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  strtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}
inline void trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  dtrtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}
inline void trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  ctrtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}
inline void trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  ztrtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}

// --- copies a triangular matrix from the standard full format (TR) to the rectangular full packed format (TF) ---
inline void trttf(char* transr, char* uplo, integer* n,  float* a, integer* lda, float* arf, integer* info)
{
  strttf_(transr, uplo, n,  a, lda, arf, info);
}
inline void trttf(char* transr, char* uplo, integer* n,  double* a, integer* lda, double* arf, integer* info)
{
  dtrttf_(transr, uplo, n,  a, lda, arf, info);
}
inline void trttf(char* transr, char* uplo, integer* n,  scomplex* a, integer* lda, scomplex* arf, integer* info)
{
  ctrttf_(transr, uplo, n,  a, lda, arf, info);
}
inline void trttf(char* transr, char* uplo, integer* n,  dcomplex* a, integer* lda, dcomplex* arf, integer* info)
{
  ztrttf_(transr, uplo, n,  a, lda, arf, info);
}

// --- copies a triangular matrix from the standard full format (TR) to the standard packed format (TP) ---
inline void trttp(char* uplo, integer* n,  float* a, integer* lda, float* ap, integer* info)
{
  strttp_(uplo, n,  a, lda, ap, info);
}
inline void trttp(char* uplo, integer* n,  double* a, integer* lda, double* ap, integer* info)
{
  dtrttp_(uplo, n,  a, lda, ap, info);
}
inline void trttp(char* uplo, integer* n,  scomplex* a, integer* lda, scomplex* ap, integer* info)
{
  ctrttp_(uplo, n,  a, lda, ap, info);
}
inline void trttp(char* uplo, integer* n,  dcomplex* a, integer* lda, dcomplex* ap, integer* info)
{
  ztrttp_(uplo, n,  a, lda, ap, info);
}

// --- reduces the M-by-N ( M<=N) real upper trapezoidal matrix A to upper triangular form---
inline void tzrzf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  stzrzf_(m, n, a, lda, tau, work, lwork, info);
}
inline void tzrzf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  dtzrzf_(m, n, a, lda, tau, work, lwork, info);
}
inline void tzrzf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  ctzrzf_(m, n, a, lda, tau, work, lwork, info);
}
inline void tzrzf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  ztzrzf_(m, n, a, lda, tau, work, lwork, info);
}

inline void tzrqf(integer* m, integer* n, float* a, integer* lda, float* tau, integer* info)
{
  printf(" Function stzrqf() has been deprecated. Please use stzrzf() instead.\n"); 
  stzrqf_(m, n, a, lda, tau, info);
}
inline void tzrqf(integer* m, integer* n, double* a, integer* lda, double* tau, integer* info)
{
  printf(" Function dtzrqf() has been deprecated. Please use dtzrzf() instead.\n"); 
  dtzrqf_(m, n, a, lda, tau, info);
}
inline void tzrqf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, integer* info)
{
  printf(" Function ctzrqf() has been deprecated. Please use ctzrzf() instead.\n"); 
  ctzrqf_(m, n, a, lda, tau, info);
}
inline void tzrqf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, integer* info)
{
  printf(" Function ztzrqf() has been deprecated. Please use ztzrzf() instead.\n"); 
  ztzrqf_(m, n, a, lda, tau, info);
}

// --- copies a triangular matrix from the standard packed format (TP) to the standard full format (TR) ---
inline void tpttr(char* uplo, integer* n,  float* ap, float* a, integer* lda, integer* info)
{
  stpttr_(uplo, n,  ap, a, lda, info);
}
inline void tpttr(char* uplo, integer* n,  double* ap, double* a, integer* lda, integer* info)
{
  dtpttr_(uplo, n,  ap, a, lda, info);
}
inline void tpttr(char* uplo, integer* n,  scomplex* ap, scomplex* a, integer* lda, integer* info)
{
  ctpttr_(uplo, n,  ap, a, lda, info);
}
inline void tpttr(char* uplo, integer* n,  dcomplex* ap, dcomplex* a, integer* lda, integer* info)
{
  ztpttr_(uplo, n,  ap, a, lda, info);
}

// --- improves the computed solution to a system of linear equations---
inline void sprfs(char* uplo, integer* n, integer* nrhs,  float* ap,  float* afp,  integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  ssprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void sprfs(char* uplo, integer* n, integer* nrhs,  double* ap,  double* afp,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  dsprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline void sprfs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  scomplex* afp,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  csprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void sprfs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* afp,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zsprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void hprfs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  scomplex* afp,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  chprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void hprfs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* afp,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zhprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}


// ---  conjugates a complex vector ---
inline void lacgv( integer* n, scomplex* x, integer*incx)
{
  clacgv_(n, x, incx); 
}
inline void lacgv( integer* n, dcomplex* x, integer*incx)
{
  zlacgv_(n, x, incx); 
}

// ---  copies all or part of a real two-dimensional array to a complex array. ---
inline void lacp2(char *uplo, integer*m, integer* n, float* a, integer*lda, scomplex* b, integer*ldb)
{
  clacp2_(uplo, m, n, a, lda,  b, ldb); 
}
inline void lacp2(char *uplo, integer*m, integer* n, double* a, integer*lda, dcomplex* b, integer*ldb)
{
  zlacp2_(uplo, m, n, a, lda,  b, ldb);
}

// ---  multiplies a complex matrix by a square real matrix ---
inline void lacrm(integer*m, integer* n, scomplex* a, integer*lda, float* b, integer*ldb, scomplex* c, integer*ldc, float* rwork)
{
  clacrm_(m, n, a, lda,  b, ldb, c, ldc, rwork); 
}
inline void lacrm(integer*m, integer* n, dcomplex* a, integer*lda, double* b, integer*ldb, dcomplex* c, integer*ldc, double* rwork)
{
  zlacrm_(m, n, a, lda,  b, ldb, c, ldc, rwork); 
}

// --- returns the value of the 1-norm, or the Frobenius norm, or the element of largest absolute value of a complex Hermitian matrix ---
inline float lanhe(char *norm, char *uplo, integer* n, scomplex* a, integer*lda, float* work)
{
  return clanhe_(norm, uplo, n, a, lda, work); 
}
inline double lanhe(char *norm, char *uplo, integer* n, dcomplex* a, integer*lda, double* work)
{
  return zlanhe_(norm, uplo, n, a, lda, work); 
}

// --- returns the value of the 1-norm, or the Frobenius norm, or the element of largest absolute value of a complex Hermitian matrix ---
inline void larcm(integer*m, integer* n, float* a, integer*lda, scomplex* b, integer*ldb, scomplex* c, integer*ldc, float* rwork)
{
  clarcm_(m, n, a, lda, b, ldb, c, ldc, rwork); 
}
inline void larcm(integer*m, integer* n, double* a, integer*lda, dcomplex* b, integer*ldb, dcomplex* c, integer*ldc, double* rwork)
{
  zlarcm_(m, n, a, lda, b, ldb, c, ldc, rwork); 
}

// --- estimates the reciprocal of the condition number of a general real matrix A ---
inline void gecon(char* norm, integer* n,  float* a, integer* lda, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  sgecon_(norm, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline void gecon(char* norm, integer* n,  double* a, integer* lda, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  dgecon_(norm, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline void gecon(char* norm, integer* n,  scomplex* a, integer* lda, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  cgecon_(norm, n,  a, lda, anorm, rcond, work, rwork, info);
}
inline void gecon(char* norm, integer* n,  dcomplex* a, integer* lda, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  zgecon_(norm, n,  a, lda, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings intended to equilibrate an M-by-N matrix A and reduce its condition number ---
inline void geequ(integer* m, integer* n,  float* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  sgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline void geequ(integer* m, integer* n,  double* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  dgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline void geequ(integer* m, integer* n,  scomplex* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  cgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline void geequ(integer* m, integer* n,  dcomplex* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  zgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}

// --- estimates the reciprocal of the condition number of a complex Hermitian matrix A---
inline void hecon(char* uplo, integer* n,  scomplex* a, integer* lda,  integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  checon_(uplo, n, a, lda,  ipiv, anorm, rcond, work, info);
}
inline void hecon(char* uplo, integer* n,  dcomplex* a, integer* lda,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  zhecon_(uplo, n, a, lda,  ipiv, anorm, rcond, work, info);
}

// --- estimates the reciprocal of the condition number---
inline void hecon_3(char* uplo, integer* n,  scomplex* a, integer* lda,  scomplex* e, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  checon_3_(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info);
}
inline void hecon_3(char* uplo, integer* n,  dcomplex* a, integer* lda,  dcomplex* e, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  zhecon_3_(uplo, n, a, lda,  e, ipiv, anorm, rcond, work, info);
}

// --- computes row and column scalings intended to equilibrate a Hermitian matrix A ---
inline void heequb(char* uplo, integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, scomplex* work, integer* info)
{
  cheequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline void heequb(char* uplo, integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, dcomplex* work, integer* info)
{
  zheequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heev_2stage(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  cheev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}
inline void heev_2stage(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zheev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heevd_2stage(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float* w, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  cheevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void heevd_2stage(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zheevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heevr_2stage(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, integer* isuppz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  cheevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void heevr_2stage(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, integer* isuppz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zheevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heevx(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  cheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline void heevx(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline void heevx_2stage(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  cheevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline void heevx_2stage(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zheevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline void hegv(integer* itype, char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  chegv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}
inline void hegv(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zhegv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline void hegv_2stage(integer* itype, char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  chegv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}
inline void hegv_2stage(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zhegv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline void hegvd(integer* itype, char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* w, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  chegvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline void hegvd(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* w, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  zhegvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes selected eigenvalues, and optionally, eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline void hegvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  chegvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline void hegvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  zhegvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- improves the computed solution to a system of linear equations ---
inline void herfs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  cherfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline void herfs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  zherfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- improves the computed solution to a system of linear equations ---
inline void herfsx(char* uplo, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  float* s,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  cherfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void herfsx(char* uplo, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  double* s,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zherfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices---
inline void hesv(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  chesv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void hesv(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zhesv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline void hesv_aa(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  chesv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void hesv_aa(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zhesv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline void hesv_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  chesv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline void hesv_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zhesv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline void hesv_rk(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  chesv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline void hesv_rk(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zhesv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline void hesvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  chesvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}
inline void hesvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  zhesvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline void hesvxx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  chesvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline void hesvxx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  zhesvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- applies an elementary permutation on the rows and columns of a Hermitian matrix ---
inline void heswapr(char* uplo, integer* n, scomplex* a, integer* lda, integer* i1, integer* i2)
{
  cheswapr_(uplo, n, a, lda, i1, i2);
}
inline void heswapr(char* uplo, integer* n, dcomplex* a, integer* lda, integer* i1, integer* i2)
{
  zheswapr_(uplo, n, a, lda, i1, i2);
}

// --- HETRD reduces a complex Hermitian matrix A to real symmetric tridiagonal form T ---
inline void hetrd(char* uplo, integer* n, scomplex* a, integer* lda, float* d, float* e, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  chetrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}
inline void hetrd(char* uplo, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  zhetrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}

// --- computes the factorization of a complex Hermitian matrix A ---
inline void hetrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  chetrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void hetrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zhetrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a complex hermitian matrix A ---
inline void hetrf_aa(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  chetrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void hetrf_aa(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zhetrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a real hermitian matrix A ---
inline void hetrf_aa_2stage(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* work, integer* lwork, integer* info)
{
  chetrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline void hetrf_aa_2stage(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* work, integer* lwork, integer* info)
{
  zhetrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}

// --- computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline void hetrf_rk(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  chetrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline void hetrf_rk(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zhetrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}

//--- computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline void hetrf_rook(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  chetrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void hetrf_rook(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zhetrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the inverse of a complex Hermitian indefinite matrix A ---
inline void hetri(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* info)
{
  chetri_(uplo, n, a, lda, ipiv, work, info);
}
inline void hetri(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* info)
{
  zhetri_(uplo, n, a, lda, ipiv, work, info);
}

// --- computes the inverse of a complex Hermitian indefinite matrix A ---
inline void hetri_3(char* uplo, integer* n, scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  chetri_3_(uplo, n, a, lda,  e,  ipiv, work, lwork, info);
}
inline void hetri_3(char* uplo, integer* n, dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zhetri_3_(uplo, n, a, lda,  e,  ipiv, work, lwork, info);
}

// --- computes the inverse of a COMPLEX hermitian indefinite matrix A ---
inline void hetri2(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  chetri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline void hetri2(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  zhetri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the inverse of a complex Hermitian indefinite matrix A ---
inline void hetri2x(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* nb, integer* info)
{
  chetri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline void hetri2x(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* nb, integer* info)
{
  zhetri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}

// --- solves a system of linear equations A*X = B with a complex Hermitian matrix A ---
inline void hetrs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  chetrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void hetrs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zhetrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A * X = B with a complex Hermitian matrix A ---
inline void hetrs_3(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  chetrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv, b, ldb, info);
}
inline void hetrs_3(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zhetrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a complex hermitian matrix A ---
inline void hetrs_aa(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  chetrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline void hetrs_aa(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  zhetrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}

// --- solves a system of linear equations A*X = B with a real hermitian matrix A ---
inline void hetrs_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, integer* info)
{
  chetrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline void hetrs_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, integer* info)
{
  zhetrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}

// --- computes the solution to a system of linear equations A * X = B for HE matrices ---
inline void hetrs_rook(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  chetrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline void hetrs_rook(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  zhetrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a complex Hermitian matrix A ---
inline void hetrs2(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* info)
{
  chetrs2_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, info);
}
inline void hetrs2(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* info)
{
  zhetrs2_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, info);
}

// --- adds two scaled sum of squares quantities ---
inline void combssq(float *v1, float *v2)
{
  scombssq_(v1, v2);
}  
inline void combssq(double *v1, double *v2)
{
  dcombssq_(v1, v2);
}

// --- forms the 1-norm of the complex vector using the true absolute value ---
inline float sum1(integer *n, scomplex *cx, integer *incx)
{
  return scsum1_(n, cx, incx);
} 
inline double sum1(integer *n, dcomplex *cx, integer *incx)
{
  return dzsum1_(n, cx, incx);
} 

// --- computes the LU factorization of a general band matrix using the unblocked version of the algorithm ---
inline void gbtf2(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, integer *ipiv, integer *info)
{
  sgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info); 
} 
inline void gbtf2(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, integer *ipiv, integer *info)
{
  dgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info);  
}
inline void gbtf2(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, integer *ipiv, integer *info)
{
  cgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info);  
}
inline void gbtf2(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, integer *ipiv, integer *info)
{
  zgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info);  
}

// --- recursively computes a LQ factorization of a float M-by-N matrix A ---
inline void gelqt3(integer *m, integer *n, float *a, integer *lda, float *t, integer *ldt, integer *info)
{
  sgelqt3_(m, n, a, lda, t, ldt, info);
}
inline void gelqt3(integer *m, integer *n, double *a, integer *lda, double *t, integer *ldt, integer *info)
{
  dgelqt3_(m, n, a, lda, t, ldt, info);
}
inline void gelqt3(integer *m, integer *n, scomplex *a, integer *lda, scomplex *t, integer *ldt, integer *info)
{
  cgelqt3_(m, n, a, lda, t, ldt, info);
}
inline void gelqt3(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, integer *info)
{
  zgelqt3_(m, n, a, lda, t, ldt, info);
}

// --- computes the QL factorization of a general rectangular matrix using an unblocked algorithm ---
inline void geql2(integer *m, integer *n, float *a, integer *lda, float *tau, float *work, integer *info)
{
  sgeql2_(m, n, a, lda, tau, work, info);
}
inline void geql2(integer *m, integer *n, double *a, integer *lda, double *tau, double *work, integer *info)
{
  dgeql2_(m, n, a, lda, tau, work, info);
}
inline void geql2(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  cgeql2_(m, n, a, lda, tau, work, info);
}
inline void geql2(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  zgeql2_(m, n, a, lda, tau, work, info);
}

// --- computes the QR factorization of a general rectangular matrix with non-negative diagonal elements using an unblocked algorithm ---
inline void geqr2p(integer *m, integer *n, float *a, integer *lda, float *tau, float *work, integer *info)
{
  sgeqr2p_(m, n, a, lda, tau, work, info);
}
inline void geqr2p(integer *m, integer *n, double *a, integer *lda, double *tau, double *work, integer *info)
{
  dgeqr2p_(m, n, a, lda, tau, work, info);
}
inline void geqr2p(integer *m, integer *n, scomplex *a, integer * lda, scomplex *tau, scomplex *work, integer *info)
{
  cgeqr2p_(m, n, a, lda, tau, work, info);
}
inline void geqr2p(integer *m, integer *n, dcomplex *a, integer * lda, dcomplex *tau, dcomplex *work, integer *info)
{
  zgeqr2p_(m, n, a, lda, tau, work, info);
}

// --- computes the RQ factorization of a general rectangular matrix using an unblocked algorithm ---
inline void gerq2(integer *m, integer *n, float *a, integer *lda, float *tau, float *work, integer *info)
{
  sgerq2_(m, n, a, lda, tau, work, info);
}
inline void gerq2(integer *m, integer *n, double *a, integer *lda, double *tau, double *work, integer *info)
{
  dgerq2_(m, n, a, lda, tau, work, info);
}
inline void gerq2(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  cgerq2_(m, n, a, lda, tau, work, info);
}
inline void gerq2(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  zgerq2_(m, n, a, lda, tau, work, info);
}

// --- solves a system of linear equations using the LU factorization with complete pivoting computed by getc2 ---
inline void gesc2(integer *n, float *a, integer *lda, float *rhs, integer *ipiv, integer *jpiv, float *scale)
{
  sgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}
inline void gesc2(integer *n, double *a, integer *lda, double *rhs, integer *ipiv, integer *jpiv, double *scale)
{
  dgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}
inline void gesc2(integer *n, scomplex *a, integer *lda, scomplex *rhs, integer *ipiv, integer *jpiv, float *scale)
{
  cgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}
inline void gesc2(integer *n, dcomplex *a, integer *lda, dcomplex *rhs, integer *ipiv, integer *jpiv, double *scale)
{
  zgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}

// --- pre-processor for the routine gesvj ---
inline void gsvj0(char *jobv, integer *m, integer *n, float *a, integer *lda, float *d, float *sva, integer *mv, float *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, float *work, integer *lwork, integer *info)
{
  sgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline void gsvj0(char *jobv, integer *m, integer *n, double *a, integer *lda, double *d, double *sva, integer *mv, double *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, double *work, integer *lwork, integer *info)
{
  dgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline void gsvj0(char *jobv, integer *m, integer *n, scomplex *a, integer *lda, scomplex *d, float *sva, integer *mv, scomplex *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, scomplex *work, integer *lwork, integer *info)
{
  cgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline void gsvj0(char *jobv, integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *d, double *sva, integer *mv, dcomplex *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, dcomplex *work, integer *lwork, integer *info)
{
  zgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}

// --- pre-processor for the routine sgesvj, applies Jacobi rotations targeting only particular pivots ---
inline void gsvj1(char *jobv, integer *m, integer *n, integer *n1, float *a, integer *lda, float *d, float *sva, integer *mv, float *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, float *work, integer *lwork, integer *info)
{
  sgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline void gsvj1(char *jobv, integer *m, integer *n, integer *n1, double *a, integer *lda, double *d, double *sva, integer *mv, double *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, double *work, integer *lwork, integer *info)
{
  dgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline void gsvj1(char *jobv, integer *m, integer *n, integer *n1, scomplex *a, integer *lda, scomplex *d, float *sva, integer *mv, scomplex *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, scomplex *work, integer *lwork, integer *info)
{
  cgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline void gsvj1(char *jobv, integer *m, integer *n, integer *n1, dcomplex *a, integer *lda, dcomplex *d, double *sva, integer *mv, dcomplex *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, dcomplex *work, integer *lwork, integer *info)
{
  zgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}

// --- solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf ---
inline void gtts2(integer *itrans, integer *n, integer *nrhs, float *dl, float *d, float *du, float *du2, integer *ipiv, float *b, integer * ldb)
{
  sgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}
inline void gtts2(integer *itrans, integer *n, integer *nrhs, double *dl, double *d, double *du, double *du2, integer *ipiv, double *b, integer * ldb)
{
  dgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}
inline void gtts2(integer *itrans, integer *n, integer *nrhs, scomplex *dl, scomplex *d, scomplex *du, scomplex *du2, integer *ipiv, scomplex *b, integer * ldb)
{
  cgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}
inline void gtts2(integer *itrans, integer *n, integer *nrhs, dcomplex *dl, dcomplex *d, dcomplex *du, dcomplex *du2, integer *ipiv, dcomplex *b, integer * ldb)
{
  zgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

// --- internal routine used by the HETRD_HB2ST subroutine ---
inline void hb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, scomplex *a, integer *lda, scomplex *v, scomplex *tau, integer *ldvt, scomplex *work)
{
  chb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
} 
inline void hb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, dcomplex *a, integer *lda, dcomplex *v, dcomplex *tau, integer *ldvt, dcomplex *work)
{
  zhb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
}

// --- estimates the reciprocal of the condition number fort HE matrices ---
inline void hecon_rook(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, float *anorm, float *rcond, scomplex *work, integer *info)
{
  checon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info); 
} 
inline void hecon_rook(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, double *anorm, double *rcond, dcomplex *work, integer *info)
{
  zhecon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info);  
}
 
// --- computes the solution to a system of linear equations A * X = B for HE matrices ---
inline void hesv_rook(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv, scomplex *b, integer *ldb, scomplex *work, integer *lwork, integer *info)
{
  chesv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline void hesv_rook(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv, dcomplex *b, integer *ldb, dcomplex *work, integer *lwork, integer *info)
{
  zhesv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the factorization of a scomplex Hermitian matrix ---
inline void hetf2(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info)
{
  chetf2_(uplo, n, a, lda, ipiv, info); 
}
inline void hetf2(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
  zhetf2_(uplo, n, a, lda, ipiv, info); 
}

// --- computes the factorization of a scomplex Hermitian indefinite matrix ---
inline void hetf2_rk(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *e, integer *ipiv, integer *info)
{
  chetf2_rk_(uplo, n, a, lda, e, ipiv, info); 
}
inline void hetf2_rk(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, integer *info)
{
  zhetf2_rk_(uplo, n, a, lda, e, ipiv, info);  
}

// --- computes the factorization of a scomplex Hermitian indefinite matrix ---
inline void hetf2_rook(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info)
{
  chetf2_rook_(uplo, n, a, lda, ipiv, info); 
}
inline void hetf2_rook(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
  zhetf2_rook_(uplo, n, a, lda, ipiv, info); 
}

// --- reduces a scomplex Hermitian matrix A to float symmetric tridiagonal form T  ---
inline void hetrd_2stage(char *vect, char *uplo, integer *n, scomplex *a, integer *lda, float *d, float *e, scomplex *tau, scomplex *hous2, integer *lhous2, scomplex *work, integer *lwork, integer *info)
{
  chetrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
}
inline void hetrd_2stage(char *vect, char *uplo, integer *n, dcomplex *a, integer *lda, double *d, double *e, dcomplex *tau, dcomplex *hous2, integer *lhous2, dcomplex *work, integer *lwork, integer *info)
{
  zhetrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
} 

// --- reduces a scomplex Hermitian band matrix A to float symmetric tridiagonal form T ---
inline void hetrd_hb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *d, float * e, scomplex *hous, integer *lhous, scomplex *work, integer *lwork, integer *info)
{
  chetrd_hb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
} 
inline void hetrd_hb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *d, double * e, dcomplex *hous, integer *lhous, dcomplex *work, integer *lwork, integer *info)
{
  zhetrd_hb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
}

// --- reduces a scomplex Hermitian matrix A to scomplex Hermitian band-diagonal form AB ---
inline void hetrd_he2hb(char *uplo, integer *n, integer *kd, scomplex *a, integer *lda, scomplex *ab, integer *ldab, scomplex *tau, scomplex *work, integer *lwork, integer *info)
{
  chetrd_he2hb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}
inline void hetrd_he2hb(char *uplo, integer *n, integer *kd, dcomplex *a, integer *lda, dcomplex *ab, integer *ldab, dcomplex *tau, dcomplex *work, integer *lwork, integer *info)
{
  zhetrd_he2hb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}

// --- computes the inverse of a scomplex Hermitian indefinite matrix A ---
inline void hetri_3x(char *uplo, integer *n, scomplex *a, integer * lda, scomplex *e, integer *ipiv, scomplex *work, integer *nb, integer * info)
{
  chetri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline void hetri_3x(char *uplo, integer *n, dcomplex *a, integer * lda, dcomplex *e, integer *ipiv, dcomplex *work, integer *nb, integer * info)
{
  zhetri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}

// --- computes the inverse of HE matrix using the factorization obtained with the bounded Bunch-Kaufman ---
inline void hetri_rook(char *uplo, integer *n, scomplex *a, integer * lda, integer *ipiv, scomplex *work, integer * info)
{
  chetri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline void hetri_rook(char *uplo, integer *n, dcomplex *a, integer * lda, integer *ipiv, dcomplex *work, integer * info)
{
  zhetri_rook_(uplo, n, a, lda, ipiv, work, info);
}

// ---  tests input for NaN ---
inline logical isnan(float *sin)
{
  return sisnan_(sin);
}
inline logical isnan(double *sin)
{
  return disnan_(sin);
}

// --- computes the Cholesky factorization of a symmetric/Hermitian positive definite band matrix (unblocked algorithm) ---
inline void pbtf2(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, integer *info)
{
  spbtf2_(uplo, n, kd, ab, ldab, info);
}
inline void pbtf2(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, integer *info)
{
  dpbtf2_(uplo, n, kd, ab, ldab, info);
}
inline void pbtf2(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, integer *info)
{
  cpbtf2_(uplo, n, kd, ab, ldab, info);
}
inline void pbtf2(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, integer *info)
{
  zpbtf2_(uplo, n, kd, ab, ldab, info);
}

// --- PSTF2 computes the Cholesky factorization with complete pivoting of a float symmetric positive semidefinite matrix ---
inline void pstf2(char *uplo, integer *n, float *a, integer *lda, integer *piv, integer *rank, float *tol, float *work, integer *info)
{
  spstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline void pstf2(char *uplo, integer *n, double *a, integer *lda, integer *piv, integer *rank, double *tol, double *work, integer *info)
{
  dpstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline void pstf2(char *uplo, integer *n, scomplex *a, integer *lda, integer *piv, integer *rank, float *tol, float *work, integer *info)
{
  cpstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline void pstf2(char *uplo, integer *n, dcomplex *a, integer *lda, integer *piv, integer *rank, double *tol, double *work, integer *info)
{
  zpstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}

// --- solves a tridiagonal system of the form AX=B using the L D LH factorization computed by pttrf ---
inline void ptts2(integer *n, integer *nrhs, float *d, float *e, float *b, integer *ldb)
{
  sptts2_(n, nrhs, d, e, b, ldb);
}
inline void ptts2(integer *n, integer *nrhs, double *d, double *e, double *b, integer *ldb)
{
  dptts2_(n, nrhs, d, e, b, ldb);
}
inline void ptts2(integer* iuplo, integer *n, integer *nrhs, float *d, scomplex *e, scomplex *b, integer *ldb)
{
  cptts2_(iuplo, n, nrhs, d, e, b, ldb);
}
inline void ptts2(integer* iuplo, integer *n, integer *nrhs, double *d, dcomplex *e, dcomplex *b, integer *ldb)
{
  zptts2_(iuplo, n, nrhs, d, e, b, ldb);
}

// --- applies a plane rotation with float cosine and scomplex sine to a pair of scomplex vectors ---
inline void rot(integer *n, scomplex *cx, integer *incx, scomplex * cy, integer *incy, float *c, scomplex *s)
{
  crot_(n, cx, incx, cy, incy, c, s);
}
inline void rot(integer *n, dcomplex *cx, integer *incx, dcomplex * cy, integer *incy, double *c, dcomplex *s)
{
  zrot_(n, cx, incx, cy, incy, c, s);
}

// --- multiplies a vector by the reciprocal of a float scalar ---
inline void rscl(integer *n, float *sa, float *sx, integer *incx)
{
  srscl_(n, sa, sx, incx);
}
inline void rscl(integer *n, double *sa, double *sx, integer *incx)
{
  drscl_(n, sa, sx, incx);
}
inline void rscl(integer *n, float *sa, scomplex *sx, integer *incx)
{
  csrscl_(n, sa, sx, incx);
}
inline void rscl(integer *n, double *sa, dcomplex *sx, integer *incx)
{
  zdrscl_(n, sa, sx, incx);
}

// --- internal routine used by the SYTRD_SB2ST subroutine ---
inline void sb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, float *a, integer *lda, float *v, float *tau, integer *ldvt, float *work)
{
  ssb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
} 
inline void sb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, double *a, integer *lda, double *v, double *tau, integer *ldvt, double *work)
{
  dsb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
}

// --- computes a matrix-vector product for scomplex vectors using a scomplex symmetric packed matrix ---
inline void spmv(char *uplo, integer *n, scomplex *alpha, scomplex *ap, scomplex *x, integer *incx, scomplex *beta, scomplex *y, integer * incy)
{
  cspmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
}
inline void spmv(char *uplo, integer *n, dcomplex *alpha, dcomplex *ap, dcomplex *x, integer *incx, dcomplex *beta, dcomplex *y, integer * incy)
{
  zspmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
}

// --- performs the symmetrical rank-1 update of a scomplex symmetric packed matrix ---
inline void spr(char *uplo, integer *n, scomplex *alpha, scomplex *x, integer *incx, scomplex *ap)
{
  cspr_(uplo, n, alpha, x, incx, ap);
}
inline void spr(char *uplo, integer *n, dcomplex *alpha, dcomplex *x, integer *incx, dcomplex *ap)
{
  zspr_(uplo, n, alpha, x, incx, ap);
}

// --- computes a matrix-vector product for a scomplex symmetric matrix ---
inline void symv(char *uplo, integer *n, scomplex *alpha, scomplex * a, integer *lda, scomplex *x, integer *incx, scomplex *beta, scomplex *y, integer *incy)
{
  csymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy); 
}
inline void symv(char *uplo, integer *n, dcomplex *alpha, dcomplex * a, integer *lda, dcomplex *x, integer *incx, dcomplex *beta, dcomplex *y, integer *incy)
{
  zsymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy); 
}

// --- reduces a float symmetric matrix A to float symmetric tridiagonal form T  ---
inline void sytrd_2stage(char *vect, char *uplo, integer *n, float *a, integer *lda, float *d, float *e, float *tau, float *hous2, integer *lhous2, float *work, integer *lwork, integer *info)
{
  ssytrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
}
inline void sytrd_2stage(char *vect, char *uplo, integer *n, double *a, integer *lda, double *d, double *e, double *tau, double *hous2, integer *lhous2, double *work, integer *lwork, integer *info)
{
  dsytrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
}

// --- reduces a float symmetric band matrix A to float symmetric tridiagonal form T ---
inline void sytrd_sb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *d, float * e, float *hous, integer *lhous, float *work, integer *lwork, integer *info)
{
  ssytrd_sb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
} 
inline void sytrd_sb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, double *ab, integer *ldab, double *d, double * e, double *hous, integer *lhous, double *work, integer *lwork, integer *info)
{
  dsytrd_sb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
}

// --- reduces a float symmetric matrix A to float symmetric band-diagonal form AB ---
inline void sytrd_sy2sb(char *uplo, integer *n, integer *kd, float *a, integer *lda, float *ab, integer *ldab, float *tau, float *work, integer *lwork, integer *info)
{
  ssytrd_sy2sb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}
inline void sytrd_sy2sb(char *uplo, integer *n, integer *kd, double *a, integer *lda, double *ab, integer *ldab, double *tau, double *work, integer *lwork, integer *info)
{
  dsytrd_sy2sb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}

// --- estimates the reciprocal of the condition number of a float symmetric matrix A ---
inline void sycon_rook(char *uplo, integer *n, float *a, integer * lda, integer *ipiv, float *anorm, float *rcond, float *work, integer * iwork, integer *info)
{
  ssycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info);
}
inline void sycon_rook(char *uplo, integer *n, double *a, integer * lda, integer *ipiv, double *anorm, double *rcond, double *work, integer * iwork, integer *info)
{
  dsycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info);
}
inline void sycon_rook(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, float *anorm, float *rcond, scomplex *work, integer *info)
{
  csycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info);
}
inline void sycon_rook(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, double *anorm, double *rcond, dcomplex *work, integer *info)
{
  zsycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info);
}

// --- converts the factorization output format ---
inline void syconvf(char *uplo, char *way, integer *n, float *a, integer *lda, float *e, integer *ipiv, integer *info)
{
  ssyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}
inline void syconvf(char *uplo, char *way, integer *n, double *a, integer *lda, double *e, integer *ipiv, integer *info)
{
  dsyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}
inline void syconvf(char *uplo, char *way, integer *n, scomplex *a, integer *lda, scomplex *e, integer *ipiv, integer *info)
{
  csyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}
inline void syconvf(char *uplo, char *way, integer *n, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, integer *info)
{
  zsyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}

// --- converts the factorization output format ---
inline void syconvf_rook(char *uplo, char *way, integer *n, float *a, integer *lda, float *e, integer *ipiv, integer *info)
{
  ssyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}
inline void syconvf_rook(char *uplo, char *way, integer *n, double *a, integer *lda, double *e, integer *ipiv, integer *info)
{
  dsyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}
inline void syconvf_rook(char *uplo, char *way, integer *n, scomplex *a, integer *lda, scomplex *e, integer *ipiv, integer *info)
{
  csyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}
inline void syconvf_rook(char *uplo, char *way, integer *n, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, integer *info)
{
  zsyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}

// --- computes the factorization of a float symmetric indefinite matrix, using the diagonal pivoting method ---
inline void sytf2(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, integer *info)
{
  ssytf2_(uplo, n, a, lda, ipiv, info);
}
inline void sytf2(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, integer *info)
{
  dsytf2_(uplo, n, a, lda, ipiv, info);
}
inline void sytf2(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info)
{
  csytf2_(uplo, n, a, lda, ipiv, info);
}
inline void sytf2(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
  zsytf2_(uplo, n, a, lda, ipiv, info);
}

// --- computes the factorization of a float symmetric indefinite matrix ---
inline void sytf2_rk(char *uplo, integer *n, float *a, integer * lda, float *e, integer *ipiv, integer *info)
{
  ssytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}
inline void sytf2_rk(char *uplo, integer *n, double *a, integer * lda, double *e, integer *ipiv, integer *info)
{
  dsytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}
inline void sytf2_rk(char *uplo, integer *n, scomplex *a, integer * lda, scomplex *e, integer *ipiv, integer *info)
{
  csytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}
inline void sytf2_rk(char *uplo, integer *n, dcomplex *a, integer * lda, dcomplex *e, integer *ipiv, integer *info)
{
  zsytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}

// --- computes the factorization of a float symmetric indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline void sytf2_rook(char *uplo, integer *n, float *a, integer * lda, integer *ipiv, integer *info)
{
  ssytf2_rook_(uplo, n, a, lda, ipiv, info);
}
inline void sytf2_rook(char *uplo, integer *n, double *a, integer * lda, integer *ipiv, integer *info)
{
  dsytf2_rook_(uplo, n, a, lda, ipiv, info);
}
inline void sytf2_rook(char *uplo, integer *n, scomplex *a, integer * lda, integer *ipiv, integer *info)
{
  csytf2_rook_(uplo, n, a, lda, ipiv, info);
}
inline void sytf2_rook(char *uplo, integer *n, dcomplex *a, integer * lda, integer *ipiv, integer *info)
{
  zsytf2_rook_(uplo, n, a, lda, ipiv, info);
}

// --- computes the inverse of a float symmetric indefinite matrix A ---
inline void sytri_3x(char *uplo, integer *n, float *a, integer * lda, float *e, integer *ipiv, float *work, integer *nb, integer *info)
{
  ssytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline void sytri_3x(char *uplo, integer *n, double *a, integer * lda, double *e, integer *ipiv, double *work, integer *nb, integer *info)
{
  dsytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline void sytri_3x(char *uplo, integer *n, scomplex *a, integer * lda, scomplex *e, integer *ipiv, scomplex *work, integer *nb, integer *info)
{
  csytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline void sytri_3x(char *uplo, integer *n, dcomplex *a, integer * lda, dcomplex *e, integer *ipiv, dcomplex *work, integer *nb, integer *info)
{
  zsytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}

// --- computes the inverse of a float symmetric matrix A ---
inline void sytri_rook(char *uplo, integer *n, float *a, integer * lda, integer *ipiv, float *work, integer *info)
{
  ssytri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline void sytri_rook(char *uplo, integer *n, double *a, integer * lda, integer *ipiv, double *work, integer *info)
{
  dsytri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline void sytri_rook(char *uplo, integer *n, scomplex *a, integer * lda, integer *ipiv, scomplex *work, integer *info)
{
  csytri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline void sytri_rook(char *uplo, integer *n, dcomplex *a, integer * lda, integer *ipiv, dcomplex *work, integer *info)
{
  zsytri_rook_(uplo, n, a, lda, ipiv, work, info);
}

// --- swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an orthogonal equivalence transformation ---
inline void tgex2(logical *wantq, logical *wantz, integer *n, float *a, integer *lda, float *b, integer *ldb, float *q, integer *ldq, float * z, integer *ldz, integer *j1, integer *n1, integer *n2, float *work, integer *lwork, integer *info)
{
  stgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info);
}
inline void tgex2(logical *wantq, logical *wantz, integer *n, double *a, integer *lda, double *b, integer *ldb, double *q, integer *ldq, double * z, integer *ldz, integer *j1, integer *n1, integer *n2, double *work, integer *lwork, integer *info)
{
  dtgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info);
}
inline void tgex2(logical *wantq, logical *wantz, integer *n, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *q, integer *ldq, scomplex *z, integer *ldz, integer *j1, integer *info)
{
  ctgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, info);
}
inline void tgex2(logical *wantq, logical *wantz, integer *n, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *q, integer *ldq, dcomplex *z, integer *ldz, integer *j1, integer *info)
{
  ztgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, info);
}

// --- solves the generalized Sylvester equation (unblocked algorithm) ---
inline void tgsy2(char *trans, integer *ijob, integer *m, integer * n, float *a, integer *lda, float *b, integer *ldb, float *c, integer * ldc, float *d, integer *ldd, float *e, integer *lde, float *f, integer *ldf, float *scale, float *rdsum, float *rdscal, integer *iwork, integer *pq, integer *info)
{
  stgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, iwork, pq, info);
}
inline void tgsy2(char *trans, integer *ijob, integer *m, integer * n, double *a, integer *lda, double *b, integer *ldb, double *c, integer * ldc, double *d, integer *ldd, double *e, integer *lde, double *f, integer *ldf, double *scale, double *rdsum, double *rdscal, integer *iwork, integer *pq, integer *info)
{
  dtgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, iwork, pq, info);
}
inline void tgsy2(char *trans, integer *ijob, integer *m, integer * n, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *c, integer *ldc, scomplex *d, integer *ldd, scomplex *e, integer *lde, scomplex *f, integer *ldf, float *scale, float *rdsum, float *rdscal, integer *info)
{
  ctgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, info);
}
inline void tgsy2(char *trans, integer *ijob, integer *m, integer * n, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *c, integer * ldc, dcomplex *d, integer *ldd, dcomplex *e, integer *lde, dcomplex *f, integer *ldf, double *scale, double *rdsum, double *rdscal, integer *info)
{
  ztgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, info);
}

// --- computes a blocked LQ factorization of a float "triangular-pentagonal" matrix C ---
inline void tplqt(integer *m, integer *n, integer *l, integer *mb, float *a, integer *lda, float *b, integer *ldb, float *t, integer *ldt, float *work, integer *info)
{
  stplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}
inline void tplqt(integer *m, integer *n, integer *l, integer *mb, double *a, integer *lda, double *b, integer *ldb, double *t, integer *ldt, double *work, integer *info)
{
  dtplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}
inline void tplqt(integer *m, integer *n, integer *l, integer *mb, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *t, integer *ldt, scomplex *work, integer *info)
{
  ctplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}
inline void tplqt(integer *m, integer *n, integer *l, integer *mb, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *t, integer *ldt, dcomplex *work, integer *info)
{
  ztplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}

// --- computes a LQ factorization of a float or scomplex "triangular-pentagonal" matrix ---
inline void tplqt2(integer *m, integer *n, integer *l, float *a, integer *lda, float *b, integer *ldb, float *t, integer *ldt, integer * info)
{
  stplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline void tplqt2(integer *m, integer *n, integer *l, double *a, integer *lda, double *b, integer *ldb, double *t, integer *ldt, integer * info)
{
  dtplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline void tplqt2(integer *m, integer *n, integer *l, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *t, integer *ldt, integer * info)
{
  ctplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline void tplqt2(integer *m, integer *n, integer *l, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *t, integer *ldt, integer * info)
{
  ztplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}

// --- applies a float orthogonal matrix Q obtained from a "triangular-pentagonal" float block reflector ---
inline void tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, float *v, integer *ldv, float *t, integer *ldt, float *a, integer *lda, float *b, integer *ldb, float * work, integer *info)
{
  stpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline void tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, double *v, integer *ldv, double *t, integer *ldt, double *a, integer *lda, double *b, integer *ldb, double * work, integer *info)
{
  dtpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline void tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, scomplex *v, integer *ldv, scomplex *t, integer *ldt, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex * work, integer *info)
{
  ctpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline void tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex * work, integer *info)
{
  ztpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}

// --- computes some or all of the right and/or left eigenvectors of a float upper quasi-triangular matrix T ---
inline void trevc3(char *side, char *howmny, logical *select, integer *n, float *t, integer *ldt, float *vl, integer *ldvl, float *vr, integer *ldvr, integer *mm, integer *m, float *work, integer *lwork, integer *info)
{
  strevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info);
}
inline void trevc3(char *side, char *howmny, logical *select, integer *n, double *t, integer *ldt, double *vl, integer *ldvl, double *vr, integer *ldvr, integer *mm, integer *m, double *work, integer *lwork, integer *info)
{
  dtrevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info);
}
inline void trevc3(char *side, char *howmny, logical *select, integer *n, scomplex *t, integer *ldt, scomplex *vl, integer *ldvl, scomplex *vr, integer *ldvr, integer *mm, integer *m, scomplex *work, integer *lwork, float *rwork, integer *lrwork, integer *info)
{
  ctrevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, rwork, lrwork, info);
}
inline void trevc3(char *side, char *howmny, logical *select, integer *n, dcomplex *t, integer *ldt, dcomplex *vl, integer *ldvl, dcomplex *vr, integer *ldvr, integer *mm, integer *m, dcomplex *work, integer *lwork, double *rwork, integer *lrwork, integer *info)
{
  ztrevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, rwork, lrwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline void orbdb1(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, integer *lwork, integer *info)
{
  sorbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void orbdb1(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  dorbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void unbdb1(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  cunbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void unbdb1(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  zunbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline void orbdb2(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, integer *lwork, integer *info)
{
  sorbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void orbdb2(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  dorbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void unbdb2(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  cunbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void unbdb2(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  zunbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline void orbdb3(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, integer *lwork, integer *info)
{
  sorbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void orbdb3(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  dorbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void unbdb3(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  cunbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline void unbdb3(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  zunbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline void orbdb4(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *phantom, float *work, integer *lwork, integer *info)
{
  sorbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}
inline void orbdb4(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *phantom, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  dorbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}
inline void unbdb4(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *phantom, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  cunbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}
inline void unbdb4(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *phantom, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  zunbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}

// --- orthogonalizes the column vector ---
inline void orbdb5(integer *m1, integer *m2, integer *n, float *x1, integer *incx1, float *x2, integer *incx2, float *q1, integer *ldq1, float *q2, integer *ldq2, float *work, integer *lwork, integer *info)
{
  sorbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline void orbdb5(integer *m1, integer *m2, integer *n, double *x1, integer *incx1, double *x2, integer *incx2, double *q1, integer *ldq1, double *q2, integer *ldq2, double *work, integer *lwork, integer *info)
{
  dorbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline void unbdb5(integer *m1, integer *m2, integer *n, scomplex *x1, integer *incx1, scomplex *x2, integer *incx2, scomplex *q1, integer *ldq1, scomplex *q2, integer *ldq2, scomplex *work, integer *lwork, integer *info)
{
  cunbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline void unbdb5(integer *m1, integer *m2, integer *n, dcomplex *x1, integer *incx1, dcomplex *x2, integer *incx2, dcomplex *q1, integer *ldq1, dcomplex *q2, integer *ldq2, dcomplex *work, integer *lwork, integer *info)
{
  zunbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}

// --- orthogonalizes the column vector ---
inline void orbdb6(integer *m1, integer *m2, integer *n, float *x1, integer *incx1, float *x2, integer *incx2, float *q1, integer *ldq1, float *q2, integer *ldq2, float *work, integer *lwork, integer *info)
{
  sorbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline void orbdb6(integer *m1, integer *m2, integer *n, double *x1, integer *incx1, double *x2, integer *incx2, double *q1, integer *ldq1, double *q2, integer *ldq2, double *work, integer *lwork, integer *info)
{
  dorbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline void unbdb6(integer *m1, integer *m2, integer *n, scomplex *x1, integer *incx1, scomplex *x2, integer *incx2, scomplex *q1, integer *ldq1, scomplex *q2, integer *ldq2, scomplex *work, integer *lwork, integer *info)
{
  cunbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline void unbdb6(integer *m1, integer *m2, integer *n, dcomplex *x1, integer *incx1, dcomplex *x2, integer *incx2, dcomplex *q1, integer *ldq1, dcomplex *q2, integer *ldq2, dcomplex *work, integer *lwork, integer *info)
{
  zunbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}

// --- generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm) ---
inline void org2l(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  sorg2l_(m, n, k, a, lda, tau, work, info);
}
inline void org2l(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  dorg2l_(m, n, k, a, lda, tau, work, info);
}
inline void ung2l(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  cung2l_(m, n, k, a, lda, tau, work, info);
}
inline void ung2l(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  zung2l_(m, n, k, a, lda, tau, work, info);
}

// --- generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf ---
inline void org2r(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  sorg2r_(m, n, k, a, lda, tau, work, info);
}
inline void org2r(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  dorg2r_(m, n, k, a, lda, tau, work, info);
}
inline void ung2r(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  cung2r_(m, n, k, a, lda, tau, work, info);
}
inline void ung2r(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  zung2r_(m, n, k, a, lda, tau, work, info);
}

// --- generates an m by n float matrix Q with orthonormal rows ---
inline void orgl2(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  sorgl2_(m, n, k, a, lda, tau, work, info);
}
inline void orgl2(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  dorgl2_(m, n, k, a, lda, tau, work, info);
}
inline void ungl2(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  cungl2_(m, n, k, a, lda, tau, work, info);
}
inline void ungl2(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  zungl2_(m, n, k, a, lda, tau, work, info);
}

// --- generates all or part of the orthogonal matrix Q from an RQ factorization determined by sgerqf ---
inline void orgr2(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  sorgr2_(m, n, k, a, lda, tau, work, info);
}
inline void orgr2(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  dorgr2_(m, n, k, a, lda, tau, work, info);
}
inline void ungr2(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  cungr2_(m, n, k, a, lda, tau, work, info);
}
inline void ungr2(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  zungr2_(m, n, k, a, lda, tau, work, info);
}

// --- generates an M-by-N float matrix Q_out with orthonormal columns ---
inline void orgtsqr(integer *m, integer *n, integer *mb, integer * nb, float *a, integer *lda, float *t, integer *ldt, float *work, integer *lwork, integer *info)
{
  sorgtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void orgtsqr(integer *m, integer *n, integer *mb, integer * nb, double *a, integer *lda, double *t, integer *ldt, double *work, integer *lwork, integer *info)
{
  dorgtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void ungtsqr(integer *m, integer *n, integer *mb, integer * nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  cungtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void ungtsqr(integer *m, integer *n, integer *mb, integer * nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  zungtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}

// --- takes an M-by-N float matrix Q_in with orthonormal columns as input and performs Householder Reconstruction ---
inline void orhr_col(integer *m, integer *n, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *d, integer *info)
{
  sorhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}
inline void orhr_col(integer *m, integer *n, integer *nb, double *a, integer *lda, double *t, integer *ldt, double *d, integer *info)
{
  dorhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}
inline void unhr_col(integer *m, integer *n, integer *nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *d, integer *info)
{
  cunhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}
inline void unhr_col(integer *m, integer *n, integer *nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *d, integer *info)
{
  zunhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a QL factorization determined by sgeqlf ---
inline void orm2l(char *side, char *trans, integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *c, integer *ldc, float *work, integer *info)
{
  sorm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void orm2l(char *side, char *trans, integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *c, integer *ldc, double *work, integer *info)
{
  dorm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void unm2l(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *c, integer *ldc, scomplex *work, integer *info)
{
  cunm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void unm2l(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work, integer *info)
{
  zunm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- multiplies a general matrix by a banded orthogonal matrix ---
inline void orm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, float *q, integer *ldq, float *c, integer * ldc, float *work, integer *lwork, integer *info)
{
  sorm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}
inline void orm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, double *q, integer *ldq, double *c, integer * ldc, double *work, integer *lwork, integer *info)
{
  dorm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}
inline void unm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, scomplex *q, integer *ldq, scomplex *c, integer * ldc, scomplex *work, integer *lwork, integer *info)
{
  cunm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}
inline void unm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, dcomplex *q, integer *ldq, dcomplex *c, integer * ldc, dcomplex *work, integer *lwork, integer *info)
{
  zunm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a RQ factorization determined by sgerqf ---
inline void ormr2(char *side, char *trans, integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *c, integer *ldc, float *work, integer *info)
{
  sormr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void ormr2(char *side, char *trans, integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *c, integer *ldc, double *work, integer *info)
{
  dormr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void unmr2(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *c, integer *ldc, scomplex *work, integer *info)
{
  cunmr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline void unmr2(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work, integer *info)
{
  zunmr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a RZ factorization determined by stzrzf ---
inline void ormr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, float *a, integer *lda, float *tau, float *c, integer *ldc, float *work, integer *info)
{
  sormr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}
inline void ormr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, double *a, integer *lda, double *tau, double *c, integer *ldc, double *work, integer *info)
{
  dormr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}
inline void unmr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, scomplex *a, integer *lda, scomplex *tau, scomplex *c, integer *ldc, scomplex *work, integer *info)
{
  cunmr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}
inline void unmr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work, integer *info)
{
  zunmr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}

// --- LA_GBAMV performs a matrix-vector operation to calculate error bounds ---
inline void la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, float *alpha, float *ab, integer *ldab, float * x, integer *incx, float *beta, float *y, integer *incy)
{
  sla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}
inline void la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha, double *ab, integer *ldab, double * x, integer *incx, double *beta, double *y, integer *incy)
{
  dla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}
inline void la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, float *alpha, scomplex *ab, integer *ldab, scomplex * x, integer *incx, float *beta, float *y, integer *incy)
{
  cla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}
inline void la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha, dcomplex *ab, integer *ldab, dcomplex * x, integer *incx, double *beta, double *y, integer *incy)
{
  zla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}

// --- estimates the Skeel condition number for a general banded matrix ---
inline float la_gbrcond(char *trans, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, float *afb, integer *ldafb, integer *ipiv, integer * cmode, float *c, integer *info, float *work, integer *iwork)
{
  return sla_gbrcond_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work, iwork);
}
inline double la_gbrcond(char *trans, integer *n, integer *kl, integer *ku, double * ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv, integer * cmode, double *c, integer *info, double *work, integer *iwork)
{
  return dla_gbrcond_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work, iwork);
}

// --- computes the infinity norm condition number of op(A)*inv(diag(c)) for general banded matrices ---
inline float la_gbrcond_c(char *trans, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float * rwork)
{
  return cla_gbrcond_c_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, capply, info, work, rwork);
}
inline double la_gbrcond_c(char *trans, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double * rwork)
{
  return zla_gbrcond_c_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, capply, info, work, rwork);
}

// --- computes the infinity norm condition number of op(A)*diag(x) for general banded matrices ---
inline float la_gbrcond_x(char *trans, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_gbrcond_x_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, rwork);
}
inline double la_gbrcond_x(char *trans, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_gbrcond_x_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, rwork);
}

// --- improves the computed solution to a system of linear equations for general banded matrices ---
inline void la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, float *ab, integer *ldab, float *afb, integer *ldafb, integer *ipiv, logical *colequ, float *c, float *b, integer *ldb, float *y, integer * ldy, float *berr_out, integer *n_norms, float *err_bnds_norm, float *err_bnds_comp, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  sla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, double *ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv, logical *colequ, double *c, double *b, integer *ldb, double *y, integer * ldy, double *berr_out, integer *n_norms, double *err_bnds_norm, double *err_bnds_comp, double *res, double *ayb, double *dy, double *y_tail, double *rcond, integer *ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  dla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer * ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex * y, integer *ldy, float *berr_out, integer *n_norms, float * err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer *ithresh, float * rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  cla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer * ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex * y, integer *ldy, double *berr_out, integer *n_norms, double * err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer *ithresh, double * rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  zla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for general matrices  ---
inline void la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, float *a, integer *lda, float * af, integer *ldaf, integer *ipiv, logical *colequ, float *c, float *b, integer *ldb, float *y, integer *ldy, float *berr_out, integer * n_norms, float *errs_n, float *errs_c, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  sla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, double *a, integer *lda, double * af, integer *ldaf, integer *ipiv, logical *colequ, double *c, double *b, integer *ldb, double *y, integer *ldy, double *berr_out, integer * n_norms, double *errs_n, double *errs_c, double *res, double *ayb, double *dy, double *y_tail, double *rcond, integer *ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  dla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer *n_norms, float *errs_n, float *errs_c, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer * ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  cla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer *n_norms, double *errs_n, double *errs_c, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer * ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  zla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for symmetric or Hermitian positive-definite matrices ---
inline void la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *af, integer * ldaf, logical *colequ, float *c, float *b, integer *ldb, float *y, integer *ldy, float *berr_out, integer *n_norms, float * err_bnds_norm, float *err_bnds_comp, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  sla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *af, integer * ldaf, logical *colequ, double *c, double *b, integer *ldb, double *y, integer *ldy, double *berr_out, integer *n_norms, double * err_bnds_norm, double *err_bnds_comp, double *res, double *ayb, double * dy, double *y_tail, double *rcond, integer *ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  dla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer *n_norms, float * err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer *ithresh, float * rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  cla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer *n_norms, double * err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer *ithresh, double * rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  zla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for symmetric indefinite matrices ---
inline void la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *af, integer * ldaf, integer *ipiv, logical *colequ, float *c, float *b, integer * ldb, float *y, integer *ldy, float *berr_out, integer *n_norms, float *err_bnds_norm, float *err_bnds_comp, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float * rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  sla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *af, integer * ldaf, integer *ipiv, logical *colequ, double *c, double *b, integer * ldb, double *y, integer *ldy, double *berr_out, integer *n_norms, double *err_bnds_norm, double *err_bnds_comp, double *res, double *ayb, double *dy, double *y_tail, double *rcond, integer *ithresh, double * rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  dla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer * n_norms, float *err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer * ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  cla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer * n_norms, double *err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer * ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  zla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for Hermitian indefinite matrices ---
inline void la_herfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer * n_norms, float *err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer * ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  cla_herfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline void la_herfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer * n_norms, double *err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer * ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  zla_herfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) for a general banded matrix ---
inline float la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, float *ab, integer *ldab, float *afb, integer *ldafb)
{
  return sla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}
inline double la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, double *ab, integer *ldab, double *afb, integer *ldafb)
{
  return dla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}
inline float la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb)
{
  return cla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}
inline double la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb)
{
  return zla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}

// --- computes a matrix-vector product using a general matrix to calculate error bounds ---
inline void la_geamv(integer *trans, integer *m, integer *n, float *alpha, float *a, integer *lda, float *x, integer *incx, float * beta, float *y, integer *incy)
{
  sla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void la_geamv(integer *trans, integer *m, integer *n, double *alpha, double *a, integer *lda, double *x, integer *incx, double * beta, double *y, integer *incy)
{
  dla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void la_geamv(integer *trans, integer *m, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x, integer *incx, float * beta, float *y, integer *incy)
{
  cla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void la_geamv(integer *trans, integer *m, integer *n, double *alpha, dcomplex *a, integer *lda, dcomplex *x, integer *incx, double * beta, double *y, integer *incy)
{
  zla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

// --- estimates the Skeel condition number for a general matrix ---
inline float la_gercond(char *trans, integer *n, float *a, integer *lda, float *af, integer *ldaf, integer *ipiv, integer *cmode, float *c, integer * info, float *work, integer *iwork)
{
  return sla_gercond_(trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}
inline  double la_gercond(char *trans, integer *n, double *a, integer *lda, double *af, integer *ldaf, integer *ipiv, integer *cmode, double *c, integer * info, double *work, integer *iwork)
{
  return dla_gercond_(trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}
inline float la_gercond_x(char *trans, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_gercond_x_(trans, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}
inline double la_gercond_x(char *trans, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_gercond_x_(trans, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}
inline float la_gercond_c(char *trans, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_gercond_c_(trans, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}
inline double la_gercond_c(char *trans, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_gercond_c_(trans, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) ---
inline float la_gerpvgrw(integer *n, integer *ncols, float *a, integer *lda, float * af, integer *ldaf)
{
  return sla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}
inline double la_gerpvgrw(integer *n, integer *ncols, double *a, integer *lda, double * af, integer *ldaf)
{
  return dla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}
inline float la_gerpvgrw(integer *n, integer *ncols, scomplex *a, integer *lda, scomplex *af, integer *ldaf)
{
  return cla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}
inline double la_gerpvgrw(integer *n, integer *ncols, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf)
{
  return zla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}

// --- computes a matrix-vector product using a Hermitian indefinite matrix to calculate error bounds ---
inline void la_heamv(integer *uplo, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x, integer *incx, float *beta, float *y, integer *incy)
{
  cla_heamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void la_heamv(integer *uplo, integer *n, double *alpha, dcomplex *a, integer *lda, dcomplex *x, integer *incx, double *beta, double *y, integer *incy)
{
  zla_heamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

// --- computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian indefinite matrices ---
inline float la_hercond_c(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_hercond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork);
}
inline double la_hercond_c(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_hercond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork);
}

// --- computes the infinity norm condition number of op(A)*diag(x) for Hermitian indefinite matrices ---
inline float la_hercond_x(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_hercond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork);
}
inline double la_hercond_x(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_hercond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork);
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) ---
inline float la_herpvgrw(char *uplo, integer *n, integer *info, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *work)
{
  return cla_herpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work); 
}
inline double la_herpvgrw(char *uplo, integer *n, integer *info, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *work)
{
  return zla_herpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work); 
}

// --- computes a component-wise relative backward error ---
inline void la_lin_berr(integer *n, integer *nz, integer *nrhs, float *res, float *ayb, float *berr)
{
  sla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}
inline void la_lin_berr(integer *n, integer *nz, integer *nrhs, double *res, double *ayb, double *berr)
{
  dla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}
inline void la_lin_berr(integer *n, integer *nz, integer *nrhs, scomplex *res, float *ayb, float *berr)
{
  cla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}
inline void la_lin_berr(integer *n, integer *nz, integer *nrhs, dcomplex *res, double *ayb, double *berr)
{
  zla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}

// --- estimates the Skeel condition number for a symmetric positive-definite matrix ---
inline float la_porcond(char *uplo, integer *n, float *a, integer *lda, float *af, integer *ldaf, integer *cmode, float *c, integer *info, float *work, integer *iwork)
{
  return sla_porcond_(uplo, n, a, lda, af, ldaf, cmode, c, info, work, iwork); 
}
inline double la_porcond(char *uplo, integer *n, double *a, integer *lda, double *af, integer *ldaf, integer *cmode, double *c, integer *info, double *work, integer *iwork)
{
  return dla_porcond_(uplo, n, a, lda, af, ldaf, cmode, c, info, work, iwork); 
}

// --- computes the infinity norm condition number of op(A)*diag(x) for Hermitian positive-definite matrices ---
inline float la_porcond_x(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_porcond_x_(uplo, n, a, lda, af, ldaf, x, info, work, rwork); 
}
inline double la_porcond_x(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_porcond_x_(uplo, n, a, lda, af, ldaf, x, info, work, rwork); 
}

// --- computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian positive-definite matrices ---
inline float la_porcond_c(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_porcond_c_(uplo, n, a, lda, af, ldaf, c, capply, info, work, rwork); 
}
inline double la_porcond_c(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_porcond_c_(uplo, n, a, lda, af, ldaf, c, capply, info, work, rwork); 
}

// --- computes a matrix-vector product using a symmetric indefinite matrix to calculate error bounds ---
inline void la_syamv(integer *uplo, integer *n, float *alpha, float *a, integer *lda, float *x, integer *incx, float *beta, float *y, integer *incy)
{
  sla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void la_syamv(integer *uplo, integer *n, double *alpha, double *a, integer *lda, double *x, integer *incx, double *beta, double *y, integer *incy)
{
  dla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void la_syamv(integer *uplo, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x, integer *incx, float *beta, float *y, integer *incy)
{
  cla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void la_syamv(integer *uplo, integer *n, double *alpha, dcomplex *a, integer *lda, dcomplex *x, integer *incx, double *beta, double *y, integer *incy)
{
  zla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Hermitian positive-definite matrix ---
inline float la_porpvgrw(char *uplo, integer *ncols, float *a, integer *lda, float * af, integer *ldaf, float *work)
{
  return sla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}
inline double la_porpvgrw(char *uplo, integer *ncols, double *a, integer *lda, double * af, integer *ldaf, double *work)
{
  return dla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}
inline float la_porpvgrw(char *uplo, integer *ncols, scomplex *a, integer *lda, scomplex * af, integer *ldaf, float *work)
{
  return cla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}
inline double la_porpvgrw(char *uplo, integer *ncols, dcomplex *a, integer *lda, dcomplex * af, integer *ldaf, double *work)
{
  return zla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}

// --- estimates the Skeel condition number for a symmetric indefinite matrix ---
inline float la_syrcond(char *uplo, integer *n, float *a, integer *lda, float *af, integer *ldaf, integer *ipiv, integer *cmode, float *c, integer * info, float *work, integer *iwork)
{
  return sla_syrcond_(uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}
inline double la_syrcond(char *uplo, integer *n, double *a, integer *lda, double *af, integer *ldaf, integer *ipiv, integer *cmode, double *c, integer * info, double *work, integer *iwork)
{
  return dla_syrcond_(uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}

// --- computes the infinity norm condition number of op(A)*diag(x) for symmetric indefinite matrices ---
inline float la_syrcond_x(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_syrcond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}
inline double la_syrcond_x(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_syrcond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}

// --- LA_SYRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for symmetric indefinite matrices ---
inline float la_syrcond_c(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_syrcond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}
inline double la_syrcond_c(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_syrcond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefinite matrix ---
inline float la_syrpvgrw(char *uplo, integer *n, integer *info, float *a, integer * lda, float *af, integer *ldaf, integer *ipiv, float *work)
{
  return sla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}
inline double la_syrpvgrw(char *uplo, integer *n, integer *info, double *a, integer * lda, double *af, integer *ldaf, integer *ipiv, double *work)
{
  return dla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}
inline float la_syrpvgrw(char *uplo, integer *n, integer *info, scomplex *a, integer * lda, scomplex *af, integer *ldaf, integer *ipiv, float *work)
{
  return cla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}
inline double la_syrpvgrw(char *uplo, integer *n, integer *info, dcomplex *a, integer * lda, dcomplex *af, integer *ldaf, integer *ipiv, double *work)
{
  return zla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}

// --- adds a vector into a doubled-single vector ---
inline void la_wwaddw(integer *n, float *x, float *y, float *w)
{
  sla_wwaddw_(n, x, y, w);
}
inline void la_wwaddw(integer *n, double *x, double *y, double *w)
{
  dla_wwaddw_(n, x, y, w);
}
inline void la_wwaddw(integer *n, scomplex *x, scomplex *y, scomplex *w)
{
  cla_wwaddw_(n, x, y, w);
}
inline void la_wwaddw(integer *n, dcomplex *x, dcomplex *y, dcomplex *w)
{
  zla_wwaddw_(n, x, y, w);
}

// --- takes as input the values computed by LAMCH for underflow and overflow, and returns the square root of each of these values ---
inline void labad(float *small_val, float *large)
{
  slabad_( small_val, large);
}
inline void labad(double *small_val, double *large)
{
  dlabad_( small_val, large);
}

// --- reduces the first nb rows and columns of a general matrix to a bidiagonal form ---
inline void labrd(integer *m, integer *n, integer *nb, float *a, integer *lda, float *d, float *e, float *tauq, float *taup, float *x, integer *ldx, float *y, integer *ldy)
{
  slabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}
inline void labrd(integer *m, integer *n, integer *nb, double *a, integer *lda, double *d, double *e, double *tauq, double *taup, double *x, integer *ldx, double *y, integer *ldy)
{
  dlabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}
inline void labrd(integer *m, integer *n, integer *nb, scomplex *a, integer *lda, float *d, float *e, scomplex *tauq, scomplex *taup, scomplex *x, integer *ldx, scomplex *y, integer *ldy)
{
  clabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}
inline void labrd(integer *m, integer *n, integer *nb, dcomplex *a, integer *lda, double *d, double *e, dcomplex *tauq, dcomplex *taup, dcomplex *x, integer *ldx, dcomplex *y, integer *ldy)
{
  zlabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}

// --- estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products ---
inline void lacon(integer *n, float *v, float *x, integer *isgn, float *est, integer *kase)
{
  slacon_(n, v, x, isgn, est, kase);
}
inline void lacon(integer *n, double *v, double *x, integer *isgn, double *est, integer *kase)
{
  dlacon_(n, v, x, isgn, est, kase);
}
inline void lacon(integer *n, scomplex *v, scomplex *x, float *est, integer *kase)
{
  clacon_(n, v, x, est, kase);
}
inline void lacon(integer *n, dcomplex *v, dcomplex *x, double *est, integer *kase)
{
  zlacon_(n, v, x, est, kase);
}

// --- performs a linear transformation of a pair of scomplex vectors ---
inline void lacrt(integer *n, scomplex *cx, integer *incx, scomplex * cy, integer *incy, scomplex *c, scomplex *s)
{
  clacrt_(n, cx, incx, cy, incy, c, s); 
}
inline void lacrt(integer *n, dcomplex *cx, integer *incx, dcomplex * cy, integer *incy, dcomplex *c, dcomplex *s)
{
  zlacrt_(n, cx, incx, cy, incy, c, s); 
}

// --- performs scomplex division in float arithmetic, avoiding unnecessary overflow ---
inline void ladiv(float *a, float *b, float *c, float *d, float *p, float *q)
{
  sladiv_(a, b, c, d, p, q);
}
inline void ladiv(double *a, double *b, double *c, double *d, double *p, double *q)
{
  dladiv_(a, b, c, d, p, q);
}

// --- computes the eigenvalues of a 2-by-2 symmetric matrix ---
inline void lae2(float *a, float *b, float *c, float *rt1, float *rt2)
{
  slae2_(a, b, c, rt1, rt2);
}
inline void lae2(double *a, double *b, double *c, double *rt1, double *rt2)
{
  dlae2_(a, b, c, rt1, rt2);
}

// --- computes the number of eigenvalues of a float symmetric tridiagonal matrix ---
inline void laebz(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp, integer *nbmin, float *abstol, float * reltol, float *pivmin, float *d, float *e, float *e2, integer *nval, float *ab, float *c, integer *mout, integer *nab, float *work, integer *iwork, integer *info)
{
  slaebz_(ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork, info);
}
inline void laebz(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp, integer *nbmin, double *abstol, double * reltol, double *pivmin, double *d, double *e, double *e2, integer *nval, double *ab, double *c, integer *mout, integer *nab, double *work, integer *iwork, integer *info)
{
  dlaebz_(ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork, info);
}

// --- computes all eigenvalues and corresponding eigenvectors of an unreduced symmetric tridiagonal matrix using the divide and conquer method ---
inline void laed0(integer *icompq, integer *qsiz, integer *n, float *d, float *e, float *q, integer *ldq, float *qstore, integer *ldqs, float *work, integer *iwork, integer *info)
{
  slaed0_(icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info);
}
inline void laed0(integer *icompq, integer *qsiz, integer *n, double *d, double *e, double *q, integer *ldq, double *qstore, integer *ldqs, double *work, integer *iwork, integer *info)
{
  dlaed0_(icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info);
}
inline void laed0(integer *qsiz, integer *n, float *d, float *e, scomplex *q, integer *ldq, scomplex *qstore, integer *ldqs, float *rwork, integer *iwork, integer *info)
{
  claed0_(qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info);
}
inline void laed0(integer *qsiz, integer *n, double *d, double *e, dcomplex *q, integer *ldq, dcomplex *qstore, integer *ldqs, double *rwork, integer *iwork, integer *info)
{
  zlaed0_(qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info);
}

// --- computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is tridiagonal ---
inline void laed1(integer *n, float *d, float *q, integer *ldq, integer *indxq, float *rho, integer *cutpnt, float *work, integer * iwork, integer *info)
{
  slaed1_(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info); 
}
inline void laed1(integer *n, double *d, double *q, integer *ldq, integer *indxq, double *rho, integer *cutpnt, double *work, integer * iwork, integer *info)
{
  dlaed1_(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info); 
}

// --- Merges eigenvalues and deflates secular equation. Used when the original matrix is tridiagonal ---
inline void laed2(integer *k, integer *n, integer *n1, float *d, float *q, integer *ldq, integer *indxq, float *rho, float *z, float * dlamda, float *w, float *q2, integer *indx, integer *indxc, integer * indxp, integer *coltyp, integer *info)
{
  slaed2_(k, n, n1, d, q, ldq, indxq, rho, z, dlamda, w, q2, indx, indxc, indxp, coltyp, info); 
}
inline void laed2(integer *k, integer *n, integer *n1, double *d, double *q, integer *ldq, integer *indxq, double *rho, double *z, double * dlamda, double *w, double *q2, integer *indx, integer *indxc, integer * indxp, integer *coltyp, integer *info)
{
  dlaed2_(k, n, n1, d, q, ldq, indxq, rho, z, dlamda, w, q2, indx, indxc, indxp, coltyp, info); 
}

// --- Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is tridiagonal ---
inline void laed3(integer *k, integer *n, integer *n1, float *d, float *q, integer *ldq, float *rho, float *dlamda, float *q2, integer * indx, integer *ctot, float *w, float *s, integer *info)
{
  slaed3_(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, info); 
}
inline void laed3(integer *k, integer *n, integer *n1, double *d, double *q, integer *ldq, double *rho, double *dlamda, double *q2, integer * indx, integer *ctot, double *w, double *s, integer *info)
{
  dlaed3_(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, info); 
}

// --- Finds a single root of the secular equation ---
inline void laed4(integer *n, integer *i, float *d, float *z, float *delta, float *rho, float *dlam, integer *info)
{
  slaed4_(n, i, d, z, delta, rho, dlam, info); 
}
inline void laed4(integer *n, integer *i, double *d, double *z, double *delta, double *rho, double *dlam, integer *info)
{
  dlaed4_(n, i, d, z, delta, rho, dlam, info); 
}

// --- Solves the 2-by-2 secular equation ---
inline void laed5(integer *i, float *d, float *z, float *delta, float *rho, float *dlam)
{
  slaed5_(i, d, z, delta, rho, dlam);
}
inline void laed5(integer *i, double *d, double *z, double *delta, double *rho, double *dlam)
{
  dlaed5_(i, d, z, delta, rho, dlam);
}

// --- computes one Newton step in solution of the secular equation ---
inline void laed6(integer *kniter, logical *orgati, float *rho, float *d, float *z, float *finit, float *tau, integer *info)
{
  slaed6_(kniter, orgati, rho, d, z, finit, tau, info);  
}
inline void laed6(integer *kniter, logical *orgati, double *rho, double *d, double *z, double *finit, double *tau, integer *info)
{
  dlaed6_(kniter, orgati, rho, d, z, finit, tau, info);  
}

// --- computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. used when the original matrix is dense ---
inline void laed7(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, float *d, float *q, integer *ldq, integer *indxq, float *rho, integer *cutpnt, float * qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, float *givnum, float *work, integer *iwork, integer *info)
{
  slaed7_(icompq, n, qsiz, tlvls, curlvl, curpbm, d, q, ldq, indxq,  rho, cutpnt, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, info); 
}
inline void laed7(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, double *d, double *q, integer *ldq, integer *indxq, double *rho, integer *cutpnt, double * qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, double *givnum, double *work, integer *iwork, integer *info)
{
  dlaed7_(icompq, n, qsiz, tlvls, curlvl, curpbm, d, q, ldq, indxq, rho, cutpnt, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, info); 
}
inline void laed7(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, float *d, scomplex *q, integer *ldq, float *rho, integer *indxq, float *qstore, integer * qptr, integer *prmptr, integer *perm, integer *givptr, integer * givcol, float *givnum, scomplex *work, float *rwork, integer *iwork, integer *info)
{
  claed7_(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info); 
}
inline void laed7(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, double *d, dcomplex *q, integer *ldq, double *rho, integer *indxq, double *qstore, integer * qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, double *givnum, dcomplex *work, double *rwork, integer *iwork, integer *info)
{
  zlaed7_(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info); 
}

// --- Merges eigenvalues and deflates secular equation. Used when the original matrix is dense ---
inline void laed8(integer *icompq, integer *k, integer *n, integer *qsiz, float *d, float *q, integer *ldq, integer *indxq, float *rho, integer *cutpnt, float *z, float *dlamda, float *q2, integer *ldq2, float *w, integer *perm, integer *givptr, integer *givcol, float * givnum, integer *indxp, integer *indx, integer *info)
{
  slaed8_(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, z, dlamda, q2, ldq2, w, perm, givptr, givcol, givnum, indxp, indx, info);
}
inline void laed8(integer *icompq, integer *k, integer *n, integer *qsiz, double *d, double *q, integer *ldq, integer *indxq, double *rho, integer *cutpnt, double *z, double *dlamda, double *q2, integer *ldq2, double *w, integer *perm, integer *givptr, integer *givcol, double * givnum, integer *indxp, integer *indx, integer *info)
{
  dlaed8_(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, z, dlamda, q2, ldq2, w, perm, givptr, givcol, givnum, indxp, indx, info);  
}
inline void laed8(integer *k, integer *n, integer *qsiz, scomplex * q, integer *ldq, float *d, float *rho, integer *cutpnt, float *z, float *dlamda, scomplex *q2, integer *ldq2, float *w, integer *indxp, integer *indx, integer *indxq, integer *perm, integer *givptr, integer *givcol, float *givnum, integer *info)
{
  claed8_(k, n, qsiz, q, ldq, d, rho, cutpnt, z, dlamda, q2, ldq2, w, indxp, indx, indxq, perm, givptr, givcol, givnum, info);
}
inline void laed8(integer *k, integer *n, integer *qsiz, dcomplex * q, integer *ldq, double *d, double *rho, integer *cutpnt, double *z, double *dlamda, dcomplex *q2, integer *ldq2, double *w, integer *indxp, integer *indx, integer *indxq, integer *perm, integer *givptr, integer *givcol, double *givnum, integer *info)
{
  zlaed8_(k, n, qsiz, q, ldq, d, rho, cutpnt, z, dlamda, q2, ldq2, w, indxp, indx, indxq, perm, givptr, givcol, givnum, info);  
}

// --- Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is dense ---
inline void laed9(integer *k, integer *kstart, integer *kstop, integer *n, float *d, float *q, integer *ldq, float *rho, float *dlamda, float *w, float *s, integer *lds, integer *info)
{
  slaed9_(k, kstart, kstop, n, d, q, ldq, rho, dlamda, w, s, lds, info); 
}
inline void laed9(integer *k, integer *kstart, integer *kstop, integer *n, double *d, double *q, integer *ldq, double *rho, double *dlamda, double *w, double *s, integer *lds, integer *info)
{
  dlaed9_(k, kstart, kstop, n, d, q, ldq, rho, dlamda, w, s, lds, info); 
}

// --- Computes the Z vector determining the rank-one modification of the diagonal matrix. Used when the original matrix is dense ---
inline void laeda(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, float *givnum, float *q, integer *qptr, float *z, float *ztemp, integer *info)
{
  slaeda_(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp, info); 
}
inline void laeda(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, double *givnum, double *q, integer *qptr, double *z, double *ztemp, integer *info)
{
  dlaeda_(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp, info); 
}

// --- computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration ---
inline void laein(logical *rightv, logical *noinit, integer *n, float *h, integer *ldh, float *wr, float *wi, float *vr, float *vi, float *b, integer *ldb, float *work, float *eps3, float *smlnum, float *bignum, integer *info)
{
  slaein_(rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum, bignum, info); 
}
inline void laein(logical *rightv, logical *noinit, integer *n, double *h, integer *ldh, double *wr, double *wi, double *vr, double *vi, double *b, integer *ldb, double *work, double *eps3, double *smlnum, double *bignum, integer *info)
{
  dlaein_(rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum, bignum, info); 
}
inline void laein(logical *rightv, logical *noinit, integer *n, scomplex *h, integer *ldh, scomplex *w, scomplex *v, scomplex *b, integer *ldb, float *rwork, float *eps3, float *smlnum, integer *info)
{
  claein_(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info); 
}
inline void laein(logical *rightv, logical *noinit, integer *n, dcomplex *h, integer *ldh, dcomplex *w, dcomplex *v, dcomplex *b, integer *ldb, double *rwork, double *eps3, double *smlnum, integer *info)
{
  zlaein_(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info); 
}

// --- computes the eigenvalues and eigenvectors of a 2-by-2 scomplex symmetric matrix ---
inline void laesy(scomplex *a, scomplex *b, scomplex *c, scomplex * rt1, scomplex *rt2, scomplex *evscal, scomplex *cs1, scomplex *sn1)
{
  claesy_(a, b, c, rt1, rt2, evscal, cs1, sn1); 
}
inline void laesy(dcomplex *a, dcomplex *b, dcomplex *c, dcomplex * rt1, dcomplex *rt2, dcomplex *evscal, dcomplex *cs1, dcomplex *sn1)
{
  zlaesy_(a, b, c, rt1, rt2, evscal, cs1, sn1); 
}

// --- LAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix ---
inline void laev2(float *a, float *b, float *c, float *rt1, float * rt2, float *cs1, float *sn1)
{
  slaev2_(a, b, c, rt1, rt2, cs1, sn1);
}
inline void laev2(double *a, double *b, double *c, double *rt1, double * rt2, double *cs1, double *sn1)
{
  dlaev2_(a, b, c, rt1, rt2, cs1, sn1);
}
inline void laev2(scomplex *a, scomplex *b, scomplex *c, float *rt1, float * rt2, float *cs1, scomplex *sn1)
{
  claev2_(a, b, c, rt1, rt2, cs1, sn1);
}
inline void laev2(dcomplex *a, dcomplex *b, dcomplex *c, double *rt1, double * rt2, double *cs1, dcomplex *sn1)
{
  zlaev2_(a, b, c, rt1, rt2, cs1, sn1);
}

// --- swaps adjacent diagonal blocks of a float upper quasi-triangular matrix in Schur canonical form, by an orthogonal similarity transformation  ---
inline void laexc(logical *wantq, integer *n, float *t, integer * ldt, float *q, integer *ldq, integer *j1, integer *n1, integer *n2, float *work, integer *info)
{
  slaexc_(wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info); 
}
inline void laexc(logical *wantq, integer *n, double *t, integer * ldt, double *q, integer *ldq, integer *j1, integer *n1, integer *n2, double *work, integer *info)
{
  dlaexc_(wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info); 
}

// --- computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, with scaling as necessary to avoid over-/underflow ---
inline void lag2(float *a, integer *lda, float *b, integer *ldb, float *safmin, float *scale1, float *scale2, float *wr1, float *wr2, float * wi)
{
  slag2_( a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi);
}
inline void lag2(double *a, integer *lda, double *b, integer *ldb, double *safmin, double *scale1, double *scale2, double *wr1, double *wr2, double * wi)
{
  dlag2_( a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi);
}

// --- computes 2-by-2 orthogonal matrices U, V, and Q, and applies them to matrices A and B such that the rows of the transformed A and B are parallel ---
inline void lags2(logical *upper, float *a1, float *a2, float *a3, float *b1, float *b2, float *b3, float *csu, float *snu, float *csv, float * snv, float *csq, float *snq)
{
  slags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}
inline void lags2(logical *upper, double *a1, double *a2, double *a3, double *b1, double *b2, double *b3, double *csu, double *snu, double *csv, double * snv, double *csq, double *snq)
{
  dlags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}
inline void lags2(logical *upper, float *a1, scomplex *a2, float *a3, float *b1, scomplex *b2, float *b3, float *csu, scomplex *snu, float *csv, scomplex *snv, float *csq, scomplex *snq)
{
  clags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}
inline void lags2(logical *upper, double *a1, dcomplex *a2, double *a3, double *b1, dcomplex *b2, double *b3, double *csu, dcomplex *snu, double *csv, dcomplex *snv, double *csq, dcomplex *snq)
{
  zlags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}

// --- computes an LU factorization of a matrix T-Î»I, where T is a general tridiagonal matrix, and Î» a scalar, using partial pivoting with row interchanges ---
inline void lagtf(integer *n, float *a, float *lambda, float *b, float *c, float *tol, float *d, integer *in, integer *info)
{
  slagtf_(n, a, lambda, b, c, tol, d, in, info);
}
inline void lagtf(integer *n, double *a, double *lambda, double *b, double *c, double *tol, double *d, integer *in, integer *info)
{
  dlagtf_(n, a, lambda, b, c, tol, d, in, info);
}

// --- performs a matrix-matrix product of the form C = Î±AB+Î²C, where A is a tridiagonal matrix, B and C are rectangular matrices, and Î± and Î² are scalars ---
inline void lagtm(char *trans, integer *n, integer *nrhs, float *alpha, float *dl, float *d, float *du, float *x, integer *ldx, float * beta, float *b, integer *ldb)
{
  return  slagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}
inline void lagtm(char *trans, integer *n, integer *nrhs, double *alpha, double *dl, double *d, double *du, double *x, integer *ldx, double *beta, double *b, integer *ldb)
{
  return  dlagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}
inline void lagtm(char *trans, integer *n, integer *nrhs, float *alpha, scomplex *dl, scomplex *d, scomplex *du, scomplex *x, integer *ldx, float *beta, scomplex *b, integer *ldb)
{
  return  clagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}
inline void lagtm(char *trans, integer *n, integer *nrhs, double *alpha, dcomplex *dl, dcomplex *d, dcomplex *du, dcomplex *x, integer *ldx, double *beta, dcomplex *b, integer *ldb)
{
  return  zlagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}

// --- solves the system of equations (T-Î»I)x = y or (T-Î»I)Tx = y ---
inline void lagts(integer *job, integer *n, float *a, float *b, float *c, float *d, integer *in, float *y, float *tol, integer *info)
{
  slagts_(job, n, a, b, c, d, in, y, tol, info);
}
inline void lagts(integer *job, integer *n, double *a, double *b, double *c, double *d, integer *in, double *y, double *tol, integer *info)
{
  dlagts_(job, n, a, b, c, d, in, y, tol, info);
}

// --- computes the Generalized Schur factorization of a float 2-by-2 matrix pencil (A,B) where B is upper triangular ---
inline void lagv2(float *a, integer *lda, float *b, integer *ldb, float *alphar, float *alphai, float *beta, float *csl, float *snl, float * csr, float *snr)
{
  slagv2_(a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr);
}
inline void lagv2(double *a, integer *lda, double *b, integer *ldb, double *alphar, double *alphai, double *beta, double *csl, double *snl, double * csr, double *snr)
{
  dlagv2_(a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr);
}

// --- computes a partial factorization of a scomplex Hermitian indefinite matrix using the Bunch-Kaufman diagonal pivoting method ---
inline void lahef(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  clahef_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline void lahef(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  zlahef_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}

// ---  factorizes a panel of a scomplex hermitian matrix A using the Aasen's algorithm ---
inline void lahef_aa(char *uplo, integer *j1, integer *m, integer *nb, scomplex *a, integer *lda, integer *ipiv, scomplex *h, integer * ldh, scomplex *work) 
{
  clahef_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work); 
}
inline void lahef_aa(char *uplo, integer *j1, integer *m, integer *nb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *h, integer * ldh, dcomplex *work) 
{
  zlahef_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work); 
}

// --- computes a partial factorization of a scomplex Hermitian indefinite matrix using bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline void lahef_rk(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, scomplex *e, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  clahef_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline void lahef_rk(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  zlahef_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}

// --- computes a partial factorization of a scomplex Hermitian matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline void lahef_rook(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  clahef_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info); 
}
inline void lahef_rook(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  zlahef_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info); 
}

// --- computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm ---
inline void lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h, integer *ldh, float *wr, float * wi, integer *iloz, integer *ihiz, float *z, integer *ldz, integer * info)
{
  slahqr_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);
}
inline void lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, double *h, integer *ldh, double *wr, double * wi, integer *iloz, integer *ihiz, double *z, integer *ldz, integer * info)
{
  dlahqr_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);
}
inline void lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, scomplex *h, integer *ldh, scomplex *w, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, integer * info)
{
  clahqr_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
}
inline void lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, dcomplex *h, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, integer * info)
{
  zlahqr_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
}

// --- reduces the specified number of first columns of a general rectangular matrix A ---
inline void lahr2(integer *n, integer *k, integer *nb, float *a, integer *lda, float *tau, float *t, integer *ldt, float *y, integer *ldy)
{
  slahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline void lahr2(integer *n, integer *k, integer *nb, double *a, integer *lda, double *tau, double *t, integer *ldt, double *y, integer *ldy)
{
  dlahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline void lahr2(integer *n, integer *k, integer *nb, scomplex *a, integer *lda, scomplex *tau, scomplex *t, integer *ldt, scomplex *y, integer *ldy)
{
  clahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline void lahr2(integer *n, integer *k, integer *nb, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *t, integer *ldt, dcomplex *y, integer *ldy)
{
  zlahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}

// --- reduces the first nb columns of a general rectangular matrix A ---
inline void lahrd(integer *n, integer *k, integer *nb, float *a, integer *lda, float *tau, float *t, integer *ldt, float *y, integer *ldy)
{
  printf(" Function slahrd() has been deprecated. Please use slahr2_() instead.\n"); 
  slahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline void lahrd(integer *n, integer *k, integer *nb, double *a, integer *lda, double *tau, double *t, integer *ldt, double *y, integer *ldy)
{
  printf(" Function dlahrd() has been deprecated. Please use dlahr2_() instead.\n"); 
  dlahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline void lahrd(integer *n, integer *k, integer *nb, scomplex *a, integer *lda, scomplex *tau, scomplex *t, integer *ldt, scomplex *y, integer *ldy)
{
  printf(" Function clahrd() has been deprecated. Please use clahr2_() instead.\n"); 
  clahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline void lahrd(integer *n, integer *k, integer *nb, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *t, integer *ldt, dcomplex *y, integer *ldy)
{
  printf(" Function zlahrd() has been deprecated. Please use zlahr2_() instead.\n"); 
  zlahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}

// --- applies one step of incremental condition estimation ---
inline void laic1(integer *job, integer *j, float *x, float *sest, float *w, float *gamma, float *sestpr, float *s, float *c__)
{
  slaic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}
inline void laic1(integer *job, integer *j, double *x, double *sest, double *w, double *gamma, double *sestpr, double *s, double *c__)
{
  dlaic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}
inline void laic1(integer *job, integer *j, scomplex *x, float *sest, scomplex *w, scomplex *gamma, float *sestpr, scomplex *s, scomplex *c__)
{
  claic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}
inline void laic1(integer *job, integer *j, dcomplex *x, double *sest, dcomplex *w, dcomplex *gamma, double *sestpr, dcomplex *s, dcomplex *c__)
{
  zlaic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}

// --- tests input for NaN by comparing two arguments for inequality  ---
inline logical laisnan(float *sin1, float *sin2)
{
  return slaisnan_(sin1, sin2);
}
inline logical laisnan(double *sin1, double *sin2)
{
  return dlaisnan_(sin1, sin2);
}

// --- solves a 1-by-1 or 2-by-2 linear system of equations of the specified form ---
inline void laln2(logical *ltrans, integer *na, integer *nw, float * smin, float *ca, float *a, integer *lda, float *d1, float *d2, float *b, integer *ldb, float *wr, float *wi, float *x, integer *ldx, float *scale, float *xnorm, integer *info)
{
  slaln2_(ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale, xnorm, info);
}
inline void laln2(logical *ltrans, integer *na, integer *nw, double * smin, double *ca, double *a, integer *lda, double *d1, double *d2, double *b, integer *ldb, double *wr, double *wi, double *x, integer *ldx, double *scale, double *xnorm, integer *info)
{
  dlaln2_(ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale, xnorm, info);
}

// ---  applies back multiplying factors in solving the least squares problem using divide and conquer SVD approach. Used by sgelsd ---
inline void lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, float *b, integer *ldb, float *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum, float *poles, float *difl, float *difr, float *z, integer *k, float *c, float *s, float *work, integer *info)
{
  slals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}
inline void lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, double *b, integer *ldb, double *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum, double *poles, double *difl, double *difr, double *z, integer *k, double *c, double *s, double *work, integer *info)
{
  dlals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}
inline void lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, scomplex *b, integer *ldb, scomplex *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum, float *poles, float * difl, float *difr, float *z, integer *k, float *c, float *s, float * work, integer *info)
{
  clals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}
inline void lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, dcomplex *b, integer *ldb, dcomplex *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum, double *poles, double * difl, double *difr, double *z, integer *k, double *c, double *s, double * work, integer *info)
{
  zlals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}

// --- computes the SVD of the coefficient matrix in compact form. Used by sgelsd ---
inline void lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, float *b, integer *ldb, float *bx, integer *ldbx, float *u, integer *ldu, float *vt, integer *k, float *difl, float *difr, float *z, float *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c, float *s, float *work, integer * iwork, integer *info)
{
  slalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}
inline void lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, double *b, integer *ldb, double *bx, integer *ldbx, double * u, integer *ldu, double *vt, integer *k, double *difl, double *difr, double * z, double *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum, double *c, double *s, double *work, integer * iwork, integer *info)
{
  dlalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}
inline void lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, scomplex *b, integer *ldb, scomplex *bx, integer *ldbx, float *u, integer *ldu, float *vt, integer *k, float *difl, float *difr, float * z, float *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c, float *s, float *work, integer * iwork, integer *info)
{
  clalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}
inline void lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, dcomplex *b, integer *ldb, dcomplex *bx, integer *ldbx, double * u, integer *ldu, double *vt, integer *k, double *difl, double *difr, double * z, double *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum, double *c, double *s, double *work, integer * iwork, integer *info)
{
  zlalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}

// --- uses the singular value decomposition of A to solve the least squares problem ---
inline void lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, float *d, float *e, float *b, integer *ldb, float *rcond, integer *rank, float *work, integer *iwork, integer *info)
{
  slalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, iwork, info);
}
inline void lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, double *d, double *e, double *b, integer *ldb, double *rcond, integer *rank, double *work, integer *iwork, integer *info)
{
  dlalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, iwork, info);
}
inline void lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, float *d, float *e, scomplex *b, integer *ldb, float *rcond, integer *rank, scomplex *work, float* rwork, integer *iwork, integer *info)
{
  clalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, rwork, iwork, info);
}
inline void lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, double *d, double *e, dcomplex *b, integer *ldb, double *rcond, integer *rank, dcomplex *work, double* rwork, integer *iwork, integer *info)
{
  zlalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, rwork, iwork, info);
}

// --- creates a permutation list to merge the entries of two independently sorted sets into a single set sorted in ascending order ---
inline void lamrg(integer *n1, integer *n2, float *a, integer * strd1, integer *strd2, integer *index)
{
  slamrg_(n1, n2, a, strd1, strd2, index);
}
inline void lamrg(integer *n1, integer *n2, double *a, integer * strd1, integer *strd2, integer *index)
{
  dlamrg_(n1, n2, a, strd1, strd2, index);
}

// --- overwrites the general float M-by-N matrix C ---
inline void lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *c, integer *ldc, float *work, integer *lwork, integer *info)
{
  slamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline void lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, double *a, integer *lda, double * t, integer *ldt, double *c, integer *ldc, double *work, integer *lwork, integer *info)
{
  dlamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline void lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, scomplex *a, integer *lda, scomplex * t, integer *ldt, scomplex *c, integer *ldc, scomplex *work, integer *lwork, integer *info)
{
  clamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline void lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, dcomplex *a, integer *lda, dcomplex * t, integer *ldt, dcomplex *c, integer *ldc, dcomplex *work, integer *lwork, integer *info)
{
  zlamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}

// --- overwrites the general float M-by-N matrix C ---
inline void lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *c, integer *ldc, float *work, integer *lwork, integer *info)
{
  slamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline void lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, double *a, integer *lda, double * t, integer *ldt, double *c, integer *ldc, double *work, integer *lwork, integer *info)
{
  dlamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline void lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, scomplex *a, integer *lda, scomplex * t, integer *ldt, scomplex *c, integer *ldc, scomplex *work, integer *lwork, integer *info)
{
  clamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline void lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, dcomplex *a, integer *lda, dcomplex * t, integer *ldt, dcomplex *c, integer *ldc, dcomplex *work, integer *lwork, integer *info)
{
  zlamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}

// --- computes the Sturm count ---
inline integer laneg(integer *n, float *d, float *lld, float *sigma, float *pivmin, integer *r__)
{
  return slaneg_(n, d, lld, sigma, pivmin, r__);
}
inline integer laneg(integer *n, double *d, double *lld, double *sigma, double *pivmin, integer *r__)
{
  return dlaneg_(n, d, lld, sigma, pivmin, r__);
}

// --- returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of general band matrix ---
inline float langb(char *norm, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, float *work)
{
  return slangb_(norm, n, kl, ku, ab, ldab, work);
}
inline double langb(char *norm, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, double *work)
{
  return dlangb_(norm, n, kl, ku, ab, ldab, work);
}
inline float langb(char *norm, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, float *work)
{
  return clangb_(norm, n, kl, ku, ab, ldab, work);
}
inline double langb(char *norm, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, double *work)
{
  return zlangb_(norm, n, kl, ku, ab, ldab, work);
}

// --- returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general tridiagonal matrix ---
inline float langt(char *norm, integer *n, float *dl, float *d, float *du)
{
  return slangt_(norm, n, dl, d, du);
}
inline double langt(char *norm, integer *n, double *dl, double *d, double *du)
{
  return dlangt_(norm, n, dl, d, du);
}
inline float langt(char *norm, integer *n, scomplex *dl, scomplex *d, scomplex *du)
{
  return clangt_(norm, n, dl, d, du);
}
inline double langt(char *norm, integer *n, dcomplex *dl, dcomplex *d, dcomplex *du)
{
  return zlangt_(norm, n, dl, d, du);
}

// --- returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a Hermitian band matrix  ---
inline float lanhb(char *norm, char *uplo, integer *n, integer *k, scomplex *ab, integer *ldab, float *work)
{
  return clanhb_(norm, uplo, n, k, ab, ldab, work); 
}
inline double lanhb(char *norm, char *uplo, integer *n, integer *k, dcomplex *ab, integer *ldab, double *work)
{
  return zlanhb_(norm, uplo, n, k, ab, ldab, work); 
}

//--- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a Hermitian matrix in RFP format ---
inline float lanhf(char *norm, char *transr, char *uplo, integer *n, scomplex *a, float *work)
{
  return clanhf_(norm, transr, uplo, n, a, work); 
}
inline double lanhf(char *norm, char *transr, char *uplo, integer *n, dcomplex *a, double *work)
{
  return zlanhf_(norm, transr, uplo, n, a, work); 
}

//--- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a scomplex Hermitian matrix supplied in packed form  ---
inline float lanhp(char *norm, char *uplo, integer *n, scomplex *ap, float *work)
{
  return clanhp_(norm, uplo, n, ap, work); 
}
inline double lanhp(char *norm, char *uplo, integer *n, dcomplex *ap, double *work)
{
  return zlanhp_(norm, uplo, n, ap, work); 
}

// --- returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix ---
inline float lanhs(char *norm, integer *n, float *a, integer *lda, float *work)
{
  return slanhs_(norm, n, a, lda, work);
}
inline double lanhs(char *norm, integer *n, double *a, integer *lda, double *work)
{
  return dlanhs_(norm, n, a, lda, work);
}
inline float lanhs(char *norm, integer *n, scomplex *a, integer *lda, float *work)
{
  return clanhs_(norm, n, a, lda, work);
}
inline double lanhs(char *norm, integer *n, dcomplex *a, integer *lda, double *work)
{
  return zlanhs_(norm, n, a, lda, work);
}

// --- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a scomplex Hermitian tridiagonal matrix  ---
inline float lanht(char *norm, integer *n, float *d, scomplex *e)
{
  return clanht_(norm, n, d, e);
}
inline double lanht(char *norm, integer *n, double *d, dcomplex *e)
{
  return zlanht_(norm, n, d, e);
}

// --- returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric band matrix  ---
inline float lansb(char *norm, char *uplo, integer *n, integer *k, float *ab, integer *ldab, float *work)
{
  return slansb_(norm, uplo, n, k, ab, ldab, work);
}
inline double lansb(char *norm, char *uplo, integer *n, integer *k, double *ab, integer *ldab, double *work)
{
  return dlansb_(norm, uplo, n, k, ab, ldab, work);
}
inline float lansb(char *norm, char *uplo, integer *n, integer *k, scomplex *ab, integer *ldab, float *work)
{
  return clansb_(norm, uplo, n, k, ab, ldab, work);
}
inline double lansb(char *norm, char *uplo, integer *n, integer *k, dcomplex *ab, integer *ldab, double *work)
{
  return zlansb_(norm, uplo, n, k, ab, ldab, work);
}

// --- returns the value of the one norm, or the infinity norm, or the element of largest absolute value of a float symmetric matrix A in RFP format  ---
inline float lansf(char *norm, char *transr, char *uplo, integer *n, float *a, float * work)
{
  return slansf_(norm, transr, uplo, n, a, work); 
}
inline double lansf(char *norm, char *transr, char *uplo, integer *n, double *a, double * work)
{
  return dlansf_(norm, transr, uplo, n, a, work); 
}

// --- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a symmetric matrix supplied in packed form ---
inline float lansp(char *norm, char *uplo, integer *n, float *ap, float *work)
{
  return slansp_(norm, uplo, n, ap, work);
}
inline double lansp(char *norm, char *uplo, integer *n, double *ap, double *work)
{
  return dlansp_(norm, uplo, n, ap, work);
}
inline float lansp(char *norm, char *uplo, integer *n, scomplex *ap, float *work)
{
  return clansp_(norm, uplo, n, ap, work);
}
inline double lansp(char *norm, char *uplo, integer *n, dcomplex *ap, double *work)
{
  return zlansp_(norm, uplo, n, ap, work);
}

// --- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a float symmetric tridiagonal matrix  ---
inline float lanst(char *norm, integer *n, float *d, float *e)
{
  return slanst_(norm, n, d, e);
}
inline double lanst(char *norm, integer *n, double *d, double *e)
{
  return dlanst_(norm, n, d, e);
}

// --- computes the Schur factorization of a float 2-by-2 nonsymmetric matrix in standard form  ---
inline void lanv2(float *a, float *b, float *c, float *d, float * rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn)
{
  slanv2_(a, b, c, d,  rt1r, rt1i, rt2r, rt2i, cs, sn);
}
inline void lanv2(double *a, double *b, double *c, double *d, double * rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn)
{
  dlanv2_(a, b, c, d,  rt1r, rt1i, rt2r, rt2i, cs, sn);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline void laorhr_col_getrfnp(integer *m, integer *n, float *a, integer *lda, float *d, integer *info)
{
  slaorhr_col_getrfnp_(m, n, a, lda, d, info);
}
inline void laorhr_col_getrfnp(integer *m, integer *n, double *a, integer *lda, double *d, integer *info)
{
  dlaorhr_col_getrfnp_(m, n, a, lda, d, info);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline void launhr_col_getrfnp(integer *m, integer *n, scomplex *a, integer *lda, scomplex *d, integer *info)
{
  claunhr_col_getrfnp_(m, n, a, lda, d, info);
}
inline void launhr_col_getrfnp(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *d, integer *info)
{
  zlaunhr_col_getrfnp_(m, n, a, lda, d, info);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline void laorhr_col_getrfnp2(integer *m, integer *n, float *a, integer *lda, float *d, integer *info)
{
  slaorhr_col_getrfnp2_(m, n, a, lda, d, info);
}
inline void laorhr_col_getrfnp2(integer *m, integer *n, double *a, integer *lda, double *d, integer *info)
{
  dlaorhr_col_getrfnp2_(m, n, a, lda, d, info);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline void launhr_col_getrfnp2(integer *m, integer *n, scomplex *a, integer *lda, scomplex *d, integer *info)
{
  claunhr_col_getrfnp2_(m, n, a, lda, d, info);
}
inline void launhr_col_getrfnp2(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *d, integer *info)
{
  zlaunhr_col_getrfnp2_(m, n, a, lda, d, info);
}

// --- measures the linear dependence of two vectors ---
inline void lapll(integer *n, float *x, integer *incx, float *y, integer *incy, float *ssmin)
{
  slapll_(n, x, incx, y, incy, ssmin);
}
inline void lapll(integer *n, double *x, integer *incx, double *y, integer *incy, double *ssmin)
{
  dlapll_(n, x, incx, y, incy, ssmin);
}
inline void lapll(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *ssmin)
{
  clapll_(n, x, incx, y, incy, ssmin);
}
inline void lapll(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *ssmin)
{
  zlapll_(n, x, incx, y, incy, ssmin);
}

// --- scales a general band matrix, using row and column scaling factors computed by sgbequ  ---
inline void laqgb(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char *equed)
{
  slaqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}
inline void laqgb(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char *equed)
{
  dlaqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}
inline void laqgb(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char *equed)
{
  claqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}
inline void laqgb(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char *equed)
{
  zlaqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}

// --- scales a general rectangular matrix, using row and column scaling factors computed by sgeequ ---
inline void laqge(integer *m, integer *n, float *a, integer *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char * equed)
{
  slaqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}
inline void laqge(integer *m, integer *n, double *a, integer *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char * equed)
{
  dlaqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}
inline void laqge(integer *m, integer *n, scomplex *a, integer *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char * equed)
{
  claqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}
inline void laqge(integer *m, integer *n, dcomplex *a, integer *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char * equed)
{
  zlaqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}

// --- scales a Hermitian band matrix, using scaling factors computed by cpbequ ---
inline void laqhb(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *s, float *scond, float *amax, char *equed)
{
  claqhb_(uplo, n, kd, ab, ldab, s, scond, amax, equed); 
}
inline void laqhb(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *s, double *scond, double *amax, char *equed)
{
  zlaqhb_(uplo, n, kd, ab, ldab, s, scond, amax, equed); 
}

// --- scales a Hermitian matrix  ---
inline void laqhe(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond, float *amax, char *equed)
{
  claqhe_(uplo, n, a, lda, s, scond, amax, equed); 
}
inline void laqhe(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond, double *amax, char *equed)
{
  zlaqhe_(uplo, n, a, lda, s, scond, amax, equed); 
}

// --- scales a Hermitian matrix stored in packed form ---
inline void laqhp(char *uplo, integer *n, scomplex *ap, float *s, float *scond, float *amax, char *equed)
{
  claqhp_(uplo, n, ap, s, scond, amax, equed); 
}
inline void laqhp(char *uplo, integer *n, dcomplex *ap, double *s, double *scond, double *amax, char *equed)
{
  zlaqhp_(uplo, n, ap, s, scond, amax, equed); 
}

// --- computes a QR factorization with column pivoting of the matrix block ---
inline void laqp2(integer *m, integer *n, integer *offset, float *a, integer *lda, integer *jpvt, float *tau, float *vn1, float *vn2, float *work)
{
  slaqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}
inline void laqp2(integer *m, integer *n, integer *offset, double *a, integer *lda, integer *jpvt, double *tau, double *vn1, double *vn2, double *work)
{
  dlaqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}
inline void laqp2(integer *m, integer *n, integer *offset, scomplex *a, integer *lda, integer *jpvt, scomplex *tau, float *vn1, float *vn2, scomplex *work)
{
  claqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}
inline void laqp2(integer *m, integer *n, integer *offset, dcomplex *a, integer *lda, integer *jpvt, dcomplex *tau, double *vn1, double *vn2, dcomplex *work)
{
  zlaqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}

// --- computes a step of QR factorization with column pivoting of a float m-by-n matrix A by using BLAS level 3 ---
inline void laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, float *a, integer *lda, integer *jpvt, float *tau, float *vn1, float *vn2, float *auxv, float *f, integer *ldf)
{
  slaqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}
inline void laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, double *a, integer *lda, integer *jpvt, double *tau, double *vn1, double *vn2, double *auxv, double *f, integer *ldf)
{
  dlaqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}
inline void laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, scomplex *a, integer *lda, integer *jpvt, scomplex *tau, float *vn1, float *vn2, scomplex *auxv, scomplex *f, integer *ldf)
{
  claqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}
inline void laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *jpvt, dcomplex *tau, double *vn1, double *vn2, dcomplex *auxv, dcomplex *f, integer *ldf)
{
  zlaqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}

// --- computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition ---
inline void laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h, integer *ldh, float *wr, float *wi, integer *iloz, integer *ihiz, float *z, integer *ldz, float *work, integer *lwork, integer *info)
{
  slaqr0_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline void laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, double *h, integer *ldh, double *wr, double *wi, integer *iloz, integer *ihiz, double *z, integer *ldz, double *work, integer *lwork, integer *info)
{
  dlaqr0_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline void laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, scomplex *h, integer *ldh, scomplex *w, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, scomplex * work, integer *lwork, integer *info)
{
  claqr0_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}
inline void laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, dcomplex *h, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, dcomplex * work, integer *lwork, integer *info)
{
  zlaqr0_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}

// --- sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts  ---
inline void laqr1(integer *n, float *h, integer *ldh, float *sr1, float *si1, float *sr2, float *si2, float *v)
{
  slaqr1_(n, h, ldh, sr1, si1, sr2, si2, v);
}
inline void laqr1(integer *n, double *h, integer *ldh, double *sr1, double *si1, double *sr2, double *si2, double *v)
{
  dlaqr1_(n, h, ldh, sr1, si1, sr2, si2, v);
}
inline void laqr1(integer *n, scomplex *h, integer *ldh, scomplex * s1, scomplex *s2, scomplex *v)
{
  claqr1_(n, h, ldh, s1, s2, v);
}
inline void laqr1(integer *n, dcomplex *h, integer *ldh, dcomplex * s1, dcomplex *s2, dcomplex *v)
{
  zlaqr1_(n, h, ldh, s1, s2, v);
}

// --- performs the orthogonal similarity transformation of a Hessenberg matrix ---
inline void laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, float *h, integer *ldh, integer *iloz, integer *ihiz, float *z, integer *ldz, integer *ns, integer *nd, float *sr, float *si, float *v, integer *ldv, integer *nh, float *t, integer *ldt, integer *nv, float *wv, integer *ldwv, float *work, integer *lwork)
{
  slaqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline void laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, double *h, integer *ldh, integer *iloz, integer *ihiz, double *z, integer *ldz, integer *ns, integer *nd, double *sr, double *si, double *v, integer *ldv, integer *nh, double *t, integer *ldt, integer *nv, double *wv, integer *ldwv, double *work, integer *lwork)
{
  dlaqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline void laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, scomplex *h, integer *ldh, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, integer *ns, integer *nd, scomplex *sh, scomplex *v, integer *ldv, integer *nh, scomplex *t, integer *ldt, integer *nv, scomplex *wv, integer *ldwv, scomplex *work, integer *lwork)
{
  claqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline void laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, dcomplex *h, integer *ldh, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, integer *ns, integer *nd, dcomplex *sh, dcomplex *v, integer *ldv, integer *nh, dcomplex *t, integer *ldt, integer *nv, dcomplex *wv, integer *ldwv, dcomplex *work, integer *lwork)
{
  zlaqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}

// --- performs the orthogonal similarity transformation of a Hessenberg matrix ---
inline void laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, float *h, integer *ldh, integer *iloz, integer *ihiz, float *z, integer *ldz, integer *ns, integer *nd, float *sr, float *si, float *v, integer *ldv, integer *nh, float *t, integer *ldt, integer *nv, float *wv, integer *ldwv, float *work, integer *lwork)
{
  slaqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline void laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, double *h, integer *ldh, integer *iloz, integer *ihiz, double *z, integer *ldz, integer *ns, integer *nd, double *sr, double *si, double *v, integer *ldv, integer *nh, double *t, integer *ldt, integer *nv, double *wv, integer *ldwv, double *work, integer *lwork)
{
  dlaqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline void laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, scomplex *h, integer *ldh, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, integer *ns, integer *nd, scomplex *sh, scomplex *v, integer *ldv, integer *nh, scomplex *t, integer *ldt, integer *nv, scomplex *wv, integer *ldwv, scomplex *work, integer *lwork)
{
  claqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline void laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, dcomplex *h, integer *ldh, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, integer *ns, integer *nd, dcomplex *sh, dcomplex *v, integer *ldv, integer *nh, dcomplex *t, integer *ldt, integer *nv, dcomplex *wv, integer *ldwv, dcomplex *work, integer *lwork)
{
  zlaqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}

// --- computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition ---
inline void laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h, integer *ldh, float *wr, float * wi, integer *iloz, integer *ihiz, float *z, integer *ldz, float *work, integer *lwork, integer *info)
{
  slaqr4_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline void laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, double *h, integer *ldh, double *wr, double * wi, integer *iloz, integer *ihiz, double *z, integer *ldz, double *work, integer *lwork, integer *info)
{
  dlaqr4_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline void laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, scomplex *h, integer *ldh, scomplex *w, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, scomplex *work, integer *lwork, integer *info)
{
  claqr4_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}
inline void laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, dcomplex *h, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, dcomplex *work, integer *lwork, integer *info)
{
  zlaqr4_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}

// --- performs a single small-bulge multi-shift QR sweep ---
inline void laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, float *sr, float *si, float *h, integer *ldh, integer *iloz, integer *ihiz, float *z, integer *ldz, float *v, integer *ldv, float *u, integer *ldu, integer *nv, float *wv, integer *ldwv, integer *nh, float *wh, integer * ldwh)
{
  slaqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}
inline void laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, double *sr, double *si, double *h, integer *ldh, integer *iloz, integer *ihiz, double *z, integer *ldz, double *v, integer *ldv, double *u, integer *ldu, integer *nv, double *wv, integer *ldwv, integer *nh, double *wh, integer * ldwh)
{
  dlaqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}
inline void laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, scomplex *s, scomplex *h, integer *ldh, integer *iloz, integer *ihiz, scomplex * z, integer *ldz, scomplex *v, integer *ldv, scomplex *u, integer *ldu, integer *nv, scomplex *wv, integer *ldwv, integer *nh, scomplex *wh, integer *ldwh)
{
  claqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}
inline void laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, dcomplex *s, dcomplex *h, integer *ldh, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, dcomplex *v, integer *ldv, dcomplex *u, integer *ldu, integer *nv, dcomplex *wv, integer *ldwv, integer *nh, dcomplex *wh, integer * ldwh)
{
  zlaqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}

// --- scales a symmetric/Hermitian band matrix, using scaling factors computed by spbequ. ---
inline void laqsb(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *s, float *scond, float *amax, char *equed)
{
  slaqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}
inline void laqsb(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, double *s, double *scond, double *amax, char *equed)
{
  dlaqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}
inline void laqsb(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *s, float *scond, float *amax, char *equed)
{
  claqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}
inline void laqsb(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *s, double *scond, double *amax, char *equed)
{
  zlaqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}

// --- scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by sppequ ---
inline void laqsp(char *uplo, integer *n, float *ap, float *s, float * scond, float *amax, char *equed)
{
  slaqsp_(uplo, n, ap, s, scond, amax, equed);
}
inline void laqsp(char *uplo, integer *n, double *ap, double *s, double * scond, double *amax, char *equed)
{
  dlaqsp_(uplo, n, ap, s, scond, amax, equed);
}
inline void laqsp(char *uplo, integer *n, scomplex *ap, float *s, float * scond, float *amax, char *equed)
{
  claqsp_(uplo, n, ap, s, scond, amax, equed);
}
inline void laqsp(char *uplo, integer *n, dcomplex *ap, double *s, double * scond, double *amax, char *equed)
{
  zlaqsp_(uplo, n, ap, s, scond, amax, equed);
}

// --- scales a symmetric/Hermitian matrix, using scaling factors computed by spoequ ---
inline void laqsy(char *uplo, integer *n, float *a, integer *lda, float *s, float *scond, float *amax, char *equed)
{
  slaqsy_(uplo, n, a, lda, s, scond, amax, equed);
}
inline void laqsy(char *uplo, integer *n, double *a, integer *lda, double *s, double *scond, double *amax, char *equed)
{
  dlaqsy_(uplo, n, a, lda, s, scond, amax, equed);
}
inline void laqsy(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond, float *amax, char *equed)
{
  claqsy_(uplo, n, a, lda, s, scond, amax, equed);
}
inline void laqsy(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond, double *amax, char *equed)
{
  zlaqsy_(uplo, n, a, lda, s, scond, amax, equed);
}

// --- solves a float quasi-triangular system of equations in float arithmetic ---
inline void laqtr(logical *ltran, logical *lreal, integer *n, float *t, integer *ldt, float *b, float *w, float *scale, float *x, float *work, integer *info)
{
  slaqtr_(ltran, lreal, n, t, ldt, b, w, scale, x, work, info);
}
inline void laqtr(logical *ltran, logical *lreal, integer *n, double *t, integer *ldt, double *b, double *w, double *scale, double *x, double *work, integer *info)
{
  dlaqtr_(ltran, lreal, n, t, ldt, b, w, scale, x, work, info);
}

// --- computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn of the tridiagonal matrix LDLT - Î»I ---
inline void lar1v(integer *n, integer *b1, integer *bn, float *lambda, float *d, float *l, float *ld, float *lld, float *pivmin, float * gaptol, float *z, logical *wantnc, integer *negcnt, float *ztz, float *mingma, integer *r, integer *isuppz, float *nrminv, float *resid, float *rqcorr, float *work)
{
  slar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}
inline void lar1v(integer *n, integer *b1, integer *bn, double *lambda, double *d, double *l, double *ld, double *lld, double *pivmin, double * gaptol, double *z, logical *wantnc, integer *negcnt, double *ztz, double * mingma, integer *r, integer *isuppz, double *nrminv, double *resid, double *rqcorr, double *work)
{
  dlar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}
inline void lar1v(integer *n, integer *b1, integer *bn, float *lambda, float *d, float *l, float *ld, float *lld, float *pivmin, float * gaptol, scomplex *z, logical *wantnc, integer *negcnt, float *ztz, float * mingma, integer *r, integer *isuppz, float *nrminv, float *resid, float *rqcorr, float *work)
{
  clar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}
inline void lar1v(integer *n, integer *b1, integer *bn, double *lambda, double *d, double *l, double *ld, double *lld, double *pivmin, double * gaptol, dcomplex *z, logical *wantnc, integer *negcnt, double *ztz, double * mingma, integer *r, integer *isuppz, double *nrminv, double *resid, double *rqcorr, double *work)
{
  zlar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}

// --- applies a vector of plane rotations with float cosines and float sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices ---
inline void lar2v(integer *n, float *x, float *y, float *z, integer *incx, float *c, float *s, integer *incc)
{
  slar2v_(n, x, y, z, incx, c, s, incc);
}
inline void lar2v(integer *n, double *x, double *y, double *z, integer *incx, double *c, double *s, integer *incc)
{
  dlar2v_(n, x, y, z, incx, c, s, incc);
}
inline void lar2v(integer *n, scomplex *x, scomplex *y, scomplex *z, integer *incx, float *c, scomplex *s, integer *incc)
{
  clar2v_(n, x, y, z, incx, c, s, incc);
}
inline void lar2v(integer *n, dcomplex *x, dcomplex *y, dcomplex *z, integer *incx, double *c, dcomplex *s, integer *incc)
{
  zlar2v_(n, x, y, z, incx, c, s, incc);
}

// --- applies an elementary reflector to a general rectangular matrix ---
inline void larf(char *side, integer *m, integer *n, float *v, integer *incv, float *tau, float *c, integer *ldc, float *work)
{
  slarf_(side, m, n, v, incv, tau, c, ldc, work);
}
inline void larf(char *side, integer *m, integer *n, double *v, integer *incv, double *tau, double *c, integer *ldc, double *work)
{
  dlarf_(side, m, n, v, incv, tau, c, ldc, work);
}
inline void larf(char *side, integer *m, integer *n, scomplex *v, integer *incv, scomplex *tau, scomplex *c, integer *ldc, scomplex *work)
{
  clarf_(side, m, n, v, incv, tau, c, ldc, work);
}
inline void larf(char *side, integer *m, integer *n, dcomplex *v, integer *incv, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work)
{
  zlarf_(side, m, n, v, incv, tau, c, ldc, work);
}

// --- applies an elementary reflector, or Householder matrix, H, to an n x n symmetric matrix C, from both the left and the right ---
inline void larfy(char *uplo, integer *n, float *v, integer *incv, float *tau, float *c, integer *ldc, float *work)
{
  slarfy_(uplo, n, v, incv, tau, c, ldc, work);
}
inline void larfy(char *uplo, integer *n, double *v, integer *incv, double *tau, double *c, integer *ldc, double *work)
{
  dlarfy_(uplo, n, v, incv, tau, c, ldc, work);
}
inline void larfy(char *uplo, integer *n, scomplex *v, integer *incv, scomplex *tau, scomplex *c, integer *ldc, scomplex *work)
{
  clarfy_(uplo, n, v, incv, tau, c, ldc, work);
}
inline void larfy(char *uplo, integer *n, dcomplex *v, integer *incv, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work)
{
  zlarfy_(uplo, n, v, incv, tau, c, ldc, work);
}

// --- generates a vector of plane rotations with float cosines and float sines ---
inline void largv(integer *n, float *x, integer *incx, float *y, integer *incy, float *c, integer *incc)
{
  slargv_(n, x, incx, y, incy, c, incc);
}
inline void largv(integer *n, double *x, integer *incx, double *y, integer *incy, double *c, integer *incc)
{
  dlargv_(n, x, incx, y, incy, c, incc);
}
inline void largv(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *c, integer *incc)
{
  clargv_(n, x, incx, y, incy, c, incc);
}
inline void largv(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *c, integer *incc)
{
  zlargv_(n, x, incx, y, incy, c, incc);
}

//--- computes the splitting points with the specified threshold ---
inline void larra(integer *n, float *d, float *e, float *e2, float * spltol, float *tnrm, integer *nsplit, integer *isplit, integer *info)
{
  slarra_(n, d, e, e2, spltol, tnrm, nsplit, isplit, info);
}
inline void larra(integer *n, double *d, double *e, double *e2, double * spltol, double *tnrm, integer *nsplit, integer *isplit, integer *info)
{
  dlarra_(n, d, e, e2, spltol, tnrm, nsplit, isplit, info);
}

//--- provides limited bisection to locate eigenvalues for more accuracy ---
inline void larrb(integer *n, float *d, float *lld, integer * ifirst, integer *ilast, float *rtol1, float *rtol2, integer *offset, float *w, float *wgap, float *werr, float *work, integer *iwork, float * pivmin, float *spdiam, integer *twist, integer *info)
{
  slarrb_(n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork, pivmin, spdiam, twist, info);
}
inline void larrb(integer *n, double *d, double *lld, integer * ifirst, integer *ilast, double *rtol1, double *rtol2, integer *offset, double *w, double *wgap, double *werr, double *work, integer *iwork, double * pivmin, double *spdiam, integer *twist, integer *info)
{
  dlarrb_(n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork, pivmin, spdiam, twist, info);
}

//--- computes the number of eigenvalues of the symmetric tridiagonal matrix ---
inline void larrc(char *jobt, integer *n, float *vl, float *vu, float *d, float *e, float *pivmin, integer *eigcnt, integer *lcnt, integer * rcnt, integer *info)
{
  slarrc_(jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info); 
}
inline void larrc(char *jobt, integer *n, double *vl, double *vu, double *d, double *e, double *pivmin, integer *eigcnt, integer *lcnt, integer * rcnt, integer *info)
{
  dlarrc_(jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info); 
}

//--- computes the eigenvalues of a symmetric tridiagonal matrix to suitable accuracy ---
inline void larrd(char *range, char *order, integer *n, float *vl, float *vu, integer *il, integer *iu, float *gers, float *reltol, float * d, float *e, float *e2, float *pivmin, integer *nsplit, integer * isplit, integer *m, float *w, float *werr, float *wl, float *wu, integer * iblock, integer *indexw, float *work, integer *iwork, integer *info)
{
  slarrd_(range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl, wu, iblock, indexw, work, iwork, info);
}
inline void larrd(char *range, char *order, integer *n, double *vl, double *vu, integer *il, integer *iu, double *gers, double *reltol, double * d, double *e, double *e2, double *pivmin, integer *nsplit, integer * isplit, integer *m, double *w, double *werr, double *wl, double *wu, integer * iblock, integer *indexw, double *work, integer *iwork, integer *info)
{
  dlarrd_(range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl, wu, iblock, indexw, work, iwork, info);
}

//--- given the tridiagonal matrix T, sets small off-diagonal elements to zero and for each unreduced block Ti, finds base representations and eigenvalues ---
inline void larre(char *range, integer *n, float *vl, float *vu, integer *il, integer *iu, float *d, float *e, float *e2, float *rtol1, float *rtol2, float *spltol, integer *nsplit, integer *isplit, integer * m, float *w, float *werr, float *wgap, integer *iblock, integer *indexw, float *gers, float *pivmin, float *work, integer *iwork, integer *info)
{
  slarre_(range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock, indexw, gers, pivmin, work, iwork, info);
}
inline void larre(char *range, integer *n, double *vl, double *vu, integer *il, integer *iu, double *d, double *e, double *e2, double *rtol1, double *rtol2, double *spltol, integer *nsplit, integer *isplit, integer * m, double *w, double *werr, double *wgap, integer *iblock, integer *indexw, double *gers, double *pivmin, double *work, integer *iwork, integer *info)
{
  dlarre_(range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock, indexw, gers, pivmin, work, iwork, info);
}

//--- finds a new relatively robust representation such that at least one of the eigenvalues is relatively isolated ---
inline void larrf(integer *n, float *d, float *l, float *ld, integer *clstrt, integer *clend, float *w, float *wgap, float *werr, float *spdiam, float *clgapl, float *clgapr, float *pivmin, float *sigma, float *dplus, float *lplus, float *work, integer *info)
{
  slarrf_(n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info);
}
inline void larrf(integer *n, double *d, double *l, double *ld, integer *clstrt, integer *clend, double *w, double *wgap, double *werr, double *spdiam, double *clgapl, double *clgapr, double *pivmin, double *sigma, double *dplus, double *lplus, double *work, integer *info)
{
  dlarrf_(n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info);
}

//--- performs refinement of the initial estimates of the eigenvalues of the matrix T ---
inline void larrj(integer *n, float *d, float *e2, integer *ifirst, integer *ilast, float *rtol, integer *offset, float *w, float *werr, float *work, integer *iwork, float *pivmin, float *spdiam, integer *info)
{
  slarrj_(n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam, info);
}
inline void larrj(integer *n, double *d, double *e2, integer *ifirst, integer *ilast, double *rtol, integer *offset, double *w, double *werr, double *work, integer *iwork, double *pivmin, double *spdiam, integer *info)
{
  dlarrj_(n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam, info);
}

//--- computes one eigenvalue of a symmetric tridiagonal matrix T to suitable accuracy ---
inline void larrk(integer *n, integer *iw, float *gl, float *gu, float *d, float *e2, float *pivmin, float *reltol, float *w, float *werr, integer *info)
{
  slarrk_(n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info);
}
inline void larrk(integer *n, integer *iw, double *gl, double *gu, double *d, double *e2, double *pivmin, double *reltol, double *w, double *werr, integer *info)
{
  dlarrk_(n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info);
}

//--- performs tests to decide whether the symmetric tridiagonal matrix T warrants expensive computations ---
inline void larrr(integer *n, float *d, float *e, integer *info)
{
  slarrr_(n, d, e, info); 
}
inline void larrr(integer *n, double *d, double *e, integer *info)
{
  dlarrr_(n, d, e, info); 
}

// --- computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenvalues of L D LT ---
inline void larrv(integer *n, float *vl, float *vu, float *d, float *l, float *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, float *wgap, integer *iblock, integer *indexw, float *gers, float *z, integer *ldz, integer *isuppz, float *work, integer *iwork, integer * info)
{
  slarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
}
inline void larrv(integer *n, double *vl, double *vu, double *d, double *l, double *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, integer *iblock, integer *indexw, double *gers, double *z, integer *ldz, integer *isuppz, double *work, integer *iwork, integer * info)
{
  dlarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
}
inline void larrv(integer *n, float *vl, float *vu, float *d, float * l, float *pivmin, integer *isplit, integer *m, integer *dol, integer * dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, float *wgap, integer *iblock, integer *indexw, float *gers, scomplex * z, integer *ldz, integer *isuppz, float *work, integer *iwork, integer *info)
{
  clarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);;
}
inline void larrv(integer *n, double *vl, double *vu, double *d, double *l, double *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, integer *iblock, integer *indexw, double *gers, dcomplex *z, integer *ldz, integer *isuppz, double *work, integer *iwork, integer * info)
{
  zlarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
}

// --- performs reciprocal diagonal scaling on a vector ---
inline void larscl2(integer *m, integer *n, float *d, float *x, integer *ldx)
{
  slarscl2_(m, n, d, x, ldx);
}
inline void larscl2(integer *m, integer *n, double *d, double *x, integer *ldx)
{
  dlarscl2_(m, n, d, x, ldx);
}
inline void larscl2(integer *m, integer *n, float *d, scomplex *x, integer *ldx)
{
  clarscl2_(m, n, d, x, ldx);
}
inline void larscl2(integer *m, integer *n, double *d, dcomplex *x, integer *ldx)
{
  zlarscl2_(m, n, d, x, ldx);
}

// --- generates a plane rotation with float cosine and float sine ---
inline void lartg(float *f, float *g, float *cs, float *sn, float *r__)
{
  slartg_( f, g, cs, sn, r__);
}
inline void lartg(double *f, double *g, double *cs, double *sn, double *r__)
{
  dlartg_( f, g, cs, sn, r__);
}
inline void lartg(scomplex *f, scomplex *g, float *cs, scomplex *sn, scomplex *r__)
{
  clartg_( f, g, cs, sn, r__);
}
inline void lartg(dcomplex *f, dcomplex *g, double *cs, dcomplex *sn, dcomplex *r__)
{
  zlartg_( f, g, cs, sn, r__);
}

// --- applies a vector of plane rotations with float cosines and float sines to the elements of a pair of vectors ---
inline void lartv(integer *n, float *x, integer *incx, float *y, integer *incy, float *c, float *s, integer *incc)
{
  slartv_(n, x, incx, y, incy, c, s, incc);
}
inline void lartv(integer *n, double *x, integer *incx, double *y, integer *incy, double *c, double *s, integer *incc)
{
  dlartv_(n, x, incx, y, incy, c, s, incc);
}
inline void lartv(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *c, scomplex *s, integer *incc)
{
  clartv_(n, x, incx, y, incy, c, s, incc);
}
inline void lartv(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *c, dcomplex *s, integer *incc)
{
  zlartv_(n, x, incx, y, incy, c, s, incc);
}

// --- returns a vector of n random float numbers from a uniform distribution ---
inline void laruv(integer *iseed, integer *n, float *x)
{
  slaruv_(iseed, n, x);  
}
inline void laruv(integer *iseed, integer *n, double *x)
{
  dlaruv_(iseed, n, x);  
}

// --- applies an elementary reflector (as returned by stzrzf) to a general matrix ---
inline void larz(char *side, integer *m, integer *n, integer *l, float *v, integer *incv, float *tau, float *c, integer *ldc, float *work)
{
  slarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}
inline void larz(char *side, integer *m, integer *n, integer *l, double *v, integer *incv, double *tau, double *c, integer *ldc, double *work)
{
  dlarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}
inline void larz(char *side, integer *m, integer *n, integer *l, scomplex *v, integer *incv, scomplex *tau, scomplex *c, integer *ldc, scomplex *work)
{
  clarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}
inline void larz(char *side, integer *m, integer *n, integer *l, dcomplex *v, integer *incv, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work)
{
  zlarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}

// --- applies a block reflector or its transpose to a general matrix ---
inline void larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, float *v, integer *ldv, float *t, integer *ldt, float *c, integer *ldc, float *work, integer *ldwork)
{
  slarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline void larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, double *v, integer *ldv, double *t, integer *ldt, double *c, integer *ldc, double *work, integer *ldwork)
{
  dlarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline void larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, scomplex *v, integer *ldv, scomplex *t, integer *ldt, scomplex *c, integer *ldc, scomplex *work, integer *ldwork)
{
  clarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline void larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt, dcomplex *c, integer *ldc, dcomplex *work, integer *ldwork)
{
  zlarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}

// --- forms the triangular factor T of a block reflector H = I - vtvH. ---
inline void larzt(char *direct, char *storev, integer *n, integer * k, float *v, integer *ldv, float *tau, float *t, integer *ldt)
{
  slarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}
inline void larzt(char *direct, char *storev, integer *n, integer * k, double *v, integer *ldv, double *tau, double *t, integer *ldt)
{
  dlarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}
inline void larzt(char *direct, char *storev, integer *n, integer * k, scomplex *v, integer *ldv, scomplex *tau, scomplex *t, integer *ldt)
{
  clarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}
inline void larzt(char *direct, char *storev, integer *n, integer * k, dcomplex *v, integer *ldv, dcomplex *tau, dcomplex *t, integer *ldt)
{
  zlarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}

// --- computes singular values of a 2-by-2 triangular matrix ---
inline void las2(float *f, float *g, float *h, float *ssmin, float * ssmax)
{
  slas2_(f, g, h, ssmin, ssmax); 
}
inline void las2(double *f, double *g, double *h, double *ssmin, double * ssmax)
{
  dlas2_(f, g, h, ssmin, ssmax); 
}

// --- performs diagonal scaling on a vector ---
inline void lascl2(integer *m, integer *n, float *d, float *x, integer *ldx)
{
  slascl2_(m, n, d, x, ldx);
}
inline void lascl2(integer *m, integer *n, double *d, double *x, integer *ldx)
{
  dlascl2_(m, n, d, x, ldx);
}
inline void lascl2(integer *m, integer *n, float *d, scomplex *x, integer *ldx)
{
  clascl2_(m, n, d, x, ldx);
}
inline void lascl2(integer *m, integer *n, double *d, dcomplex *x, integer *ldx)
{
  zlascl2_(m, n, d, x, ldx);
}

// --- computes the singular values of a float upper bidiagonal n-by-m matrix B with diagonal d and off-diagonal e ---
inline void lasd0(integer *n, integer *sqre, float *d, float *e, float *u, integer *ldu, float *vt, integer *ldvt, integer *smlsiz, integer *iwork, float *work, integer *info)
{
  slasd0_(n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info);
}
inline void lasd0(integer *n, integer *sqre, double *d, double *e, double *u, integer *ldu, double *vt, integer *ldvt, integer *smlsiz, integer *iwork, double *work, integer *info)
{
  dlasd0_(n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info); 
}

// --- computes the SVD of an upper bidiagonal matrix B of the specified size ---
inline void lasd1(integer *nl, integer *nr, integer *sqre, float * d, float *alpha, float *beta, float *u, integer *ldu, float *vt, integer *ldvt, integer *idxq, integer *iwork, float *work, integer * info)
{
  slasd1_(nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt, idxq, iwork, work, info);
}
inline void lasd1(integer *nl, integer *nr, integer *sqre, double * d, double *alpha, double *beta, double *u, integer *ldu, double *vt, integer *ldvt, integer *idxq, integer *iwork, double *work, integer * info)
{
  dlasd1_(nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt, idxq, iwork, work, info); 
}

// --- merges the two sets of singular values together into a single sorted set ---
inline void lasd2(integer *nl, integer *nr, integer *sqre, integer *k, float *d, float *z, float *alpha, float *beta, float *u, integer * ldu, float *vt, integer *ldvt, float *dsigma, float *u2, integer *ldu2, float *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info)
{
  slasd2_(nl, nr, sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2, ldvt2, idxp, idx, idxc, idxq, coltyp, info); 
}
inline void lasd2(integer *nl, integer *nr, integer *sqre, integer *k, double *d, double *z, double *alpha, double *beta, double *u, integer * ldu, double *vt, integer *ldvt, double *dsigma, double *u2, integer *ldu2, double *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info)
{
  dlasd2_(nl, nr, sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2, ldvt2, idxp, idx, idxc, idxq, coltyp, info); 
}

// --- finds all square roots of the roots of the secular equation ---
inline void lasd3(integer *nl, integer *nr, integer *sqre, integer *k, float *d, float *q, integer *ldq, float *dsigma, float *u, integer * ldu, float *u2, integer *ldu2, float *vt, integer *ldvt, float *vt2, integer *ldvt2, integer *idxc, integer *ctot, float *z, integer * info)
{
  slasd3_(nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2, ldu2, vt, ldvt, vt2, ldvt2, idxc, ctot, z, info);
}
inline void lasd3(integer *nl, integer *nr, integer *sqre, integer *k, double *d, double *q, integer *ldq, double *dsigma, double *u, integer * ldu, double *u2, integer *ldu2, double *vt, integer *ldvt, double *vt2, integer *ldvt2, integer *idxc, integer *ctot, double *z, integer * info)
{
  dlasd3_(nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2, ldu2, vt, ldvt, vt2, ldvt2, idxc, ctot, z, info); 
}

// --- computes the square root of the i-th updated eigenvalue of a positive symmetric rank-one modification to a positive diagonal matrix ---
inline void lasd4(integer *n, integer *i, float *d, float *z, float *delta, float *rho, float *sigma, float *work, integer *info)
{
  slasd4_(n, i, d, z, delta, rho, sigma, work, info);
}
inline void lasd4(integer *n, integer *i, double *d, double *z, double *delta, double *rho, double *sigma, double *work, integer *info)
{
  dlasd4_(n, i, d, z, delta, rho, sigma, work, info); 
}

// --- computes the square root of the i-th eigenvalue of a positive symmetric rank-one modification of a 2-by-2 diagonal matrix ---
inline void lasd5(integer *i, float *d, float *z, float *delta, float *rho, float *dsigma, float *work)
{
  slasd5_(i, d, z, delta, rho, dsigma, work); 
}
inline void lasd5(integer *i, double *d, double *z, double *delta, double *rho, double *dsigma, double *work)
{
  dlasd5_(i, d, z, delta, rho, dsigma, work);  
}

// --- computes the SVD of an updated upper bidiagonal matrix obtained by merging two smaller ones by appending a row ---
inline void lasd6(integer *icompq, integer *nl, integer *nr, integer *sqre, float *d, float *vf, float *vl, float *alpha, float *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum, float *poles, float * difl, float *difr, float *z, integer *k, float *c, float *s, float * work, integer *iwork, integer *info)
{
  slasd6_(icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, iwork, info); 
}
inline void lasd6(integer *icompq, integer *nl, integer *nr, integer *sqre, double *d, double *vf, double *vl, double *alpha, double *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum, double *poles, double * difl, double *difr, double *z, integer *k, double *c, double *s, double * work, integer *iwork, integer *info)
{
  dlasd6_(icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, iwork, info); 
}

// --- merges the two sets of singular values together into a single sorted set. Then it tries to deflate the size of the problem ---
inline void lasd7(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, float *d, float *z, float *zw, float *vf, float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, integer * givptr, integer *givcol, integer *ldgcol, float *givnum, integer * ldgnum, float *c, float *s, integer *info)
{
  slasd7_(icompq, nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma, idx, idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info); 
}
inline void lasd7(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, double *d, double *z, double *zw, double *vf, double *vfw, double *vl, double *vlw, double *alpha, double *beta, double *dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, integer * givptr, integer *givcol, integer *ldgcol, double *givnum, integer * ldgnum, double *c, double *s, integer *info)
{
  dlasd7_(icompq, nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma, idx, idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info); 
}

// --- finds the square roots of the roots of the secular equation, and stores, for each element in D, the distance to its two nearest poles ---
inline void lasd8(integer *icompq, integer *k, float *d, float * z, float *vf, float *vl, float *difl, float *difr, integer *lddifr, float *dsigma, float *work, integer *info)
{
  slasd8_(icompq, k, d, z, vf, vl, difl, difr, lddifr, dsigma, work, info); 
}
inline void lasd8(integer *icompq, integer *k, double *d, double * z, double *vf, double *vl, double *difl, double *difr, integer *lddifr, double *dsigma, double *work, integer *info)
{
  dlasd8_(icompq, k, d, z, vf, vl, difl, difr, lddifr, dsigma, work, info); 
}

// --- computes the singular value decomposition (SVD) of a float upper bidiagonal matrix with diagonal d and off-diagonal e ---
inline void lasda(integer *icompq, integer *smlsiz, integer *n, integer *sqre, float *d, float *e, float *u, integer *ldu, float *vt, integer *k, float *difl, float *difr, float *z, float *poles, integer * givptr, integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c, float *s, float *work, integer *iwork, integer *info)
{
  slasda_(icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info); 
}
inline void lasda(integer *icompq, integer *smlsiz, integer *n, integer *sqre, double *d, double *e, double *u, integer *ldu, double *vt, integer *k, double *difl, double *difr, double *z, double *poles, integer * givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum, double *c, double *s, double *work, integer *iwork, integer *info)
{
  dlasda_(icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info); 
}

// --- computes the SVD of a float bidiagonal matrix with diagonal d and off-diagonal e ---
inline void lasdq(char *uplo, integer *sqre, integer *n, integer * ncvt, integer *nru, integer *ncc, float *d, float *e, float *vt, integer *ldvt, float *u, integer *ldu, float *c, integer *ldc, float * work, integer *info)
{
  slasdq_(uplo, sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info); 
}
inline void lasdq(char *uplo, integer *sqre, integer *n, integer * ncvt, integer *nru, integer *ncc, double *d, double *e, double *vt, integer *ldvt, double *u, integer *ldu, double *c, integer *ldc, double * work, integer *info)
{
  dlasdq_(uplo, sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info); 
}

// --- computes the singular values of a float square bidiagonal matrix ---
inline void lasq1(integer *n, float *d, float *e, float *work, integer *info)
{
  slasq1_(n, d, e, work, info);
}
inline void lasq1(integer *n, double *d, double *e, double *work, integer *info)
{
  dlasq1_(n, d, e, work, info);
}

// --- computes all the eigenvalues of the symmetric positive definite tridiagonal matrix associated with the qd Array Z to high relative accuracy. ---
inline void lasq2(integer *n, float *z, integer *info)
{
  slasq2_(n, z, info);
}
inline void lasq2(integer *n, double *z, integer *info)
{
  dlasq2_(n, z, info);
}

// --- checks for deflation, computes a shift and calls dqds ---
inline void lasq3(integer *i0, integer *n0, float *z, integer *pp, float *dmin, float *sigma, float *desig, float *qmax, integer *nfail, integer *iter, integer *ndiv, logical *ieee, integer *ttype, float * dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *g, float * tau)
{
  slasq3_(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau);
}
inline void lasq3(integer *i0, integer *n0, double *z, integer *pp, double *dmin, double *sigma, double *desig, double *qmax, integer *nfail, integer *iter, integer *ndiv, logical *ieee, integer *ttype, double * dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *g, double * tau)
{
  dlasq3_(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau);
}

// --- computes an approximation to the smallest eigenvalue using values of d from the previous transform  ---
inline void lasq4(integer *i0, integer *n0, float *z, integer *pp, integer *n0in, float *dmin, float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau, integer *ttype, float *g)
{
  slasq4_(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
}
inline void lasq4(integer *i0, integer *n0, double *z, integer *pp, integer *n0in, double *dmin, double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau, integer *ttype, double *g)
{
  dlasq4_(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
}

// --- computes one dqds transform in ping-pong form ---
inline void lasq5(integer *i0, integer *n0, float *z, integer *pp, float *tau, float *sigma, float *dmin, float *dmin1, float *dmin2, float *dn, float *dnm1, float *dnm2, logical *ieee, float *eps)
{
  slasq5_(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps);
}
inline void lasq5(integer *i0, integer *n0, double *z, integer *pp, double *tau, double *sigma, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2, logical *ieee, double *eps)
{
  dlasq5_(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps);
}

// --- computes one dqd transform in ping-pong form ---
inline void lasq6(integer *i0, integer *n0, float *z, integer *pp, float *dmin, float *dmin1, float *dmin2, float *dn, float *dnm1, float * dnm2) 
{
  slasq6_(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2);
}
inline void lasq6(integer *i0, integer *n0, double *z, integer *pp, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double * dnm2)
{
  dlasq6_(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2);
}

// --- applies a sequence of plane rotations to a general rectangular matrix ---
inline void lasr(char *side, char *pivot, char *direct, integer *m, integer *n, float *c, float *s, float *a, integer *lda)
{
  slasr_(side, pivot, direct, m, n, c, s, a, lda);
}
inline void lasr_(char *side, char *pivot, char *direct, integer *m, integer *n, double *c, double *s, double *a, integer *lda)
{
  dlasr_(side, pivot, direct, m, n, c, s, a, lda);
}
inline void lasr(char *side, char *pivot, char *direct, integer *m, integer *n, float *c, float *s, scomplex *a, integer *lda)
{
  clasr_(side, pivot, direct, m, n, c, s, a, lda);
}
inline void lasr(char *side, char *pivot, char *direct, integer *m, integer *n, double *c, double *s, dcomplex *a, integer *lda)
{
  zlasr_(side, pivot, direct, m, n, c, s, a, lda);
}

// --- computes the singular value decomposition of a 2-by-2 triangular matrix ---
inline void lasv2(float *f, float *g, float *h, float *ssmin, float * ssmax, float *snr, float *csr, float *snl, float *csl)
{
  slasv2_( f, g, h, ssmin, ssmax, snr, csr, snl, csl);
}
inline void lasv2(double *f, double *g, double *h, double *ssmin, double * ssmax, double *snr, double *csr, double *snl, double *csl)
{
  dlasv2_( f, g, h, ssmin, ssmax, snr, csr, snl, csl);
}

// --- computes a blocked Tall-Skinny LQ factorization of a float M-by-N matrix A for M <= N ---
inline void laswlq(integer *m, integer *n, integer *mb, integer * nb, float *a, integer *lda, float *t, integer *ldt, float *work, integer *lwork, integer *info)
{
  slaswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void laswlq(integer *m, integer *n, integer *mb, integer * nb, double *a, integer *lda, double *t, integer *ldt, double *work, integer *lwork, integer *info)
{
  dlaswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void laswlq(integer *m, integer *n, integer *mb, integer * nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  claswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void laswlq(integer *m, integer *n, integer *mb, integer * nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  zlaswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}

// --- solves the Sylvester matrix equation where the matrices are of order 1 or 2 ---
inline void lasy2(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2, float *tl, integer *ldtl, float *tr, integer * ldtr, float *b, integer *ldb, float *scale, float *x, integer *ldx, float *xnorm, integer *info)
{
  slasy2_(ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info);
}
inline void lasy2(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2, double *tl, integer *ldtl, double *tr, integer * ldtr, double *b, integer *ldb, double *scale, double *x, integer *ldx, double *xnorm, integer *info)
{
  dlasy2_(ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info);
}

// --- computes a partial factorization of a float symmetric matrix using the Bunch-Kaufman diagonal pivoting method ---
inline void lasyf(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda, integer *ipiv, float *w, integer *ldw, integer *info)
{
  slasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline void lasyf(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda, integer *ipiv, double *w, integer *ldw, integer *info)
{
  dlasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline void lasyf(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  clasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline void lasyf(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  zlasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}

// --- factorizes a panel of a float symmetric matrix A using the Aasen's algorithm ---
inline void lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, float *a, integer *lda, integer *ipiv, float *h, integer *ldh, float *work)
{
  slasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}
inline void lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, double *a, integer *lda, integer *ipiv, double *h, integer *ldh, double *work)
{
  dlasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}
inline void lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, scomplex *a, integer *lda, integer *ipiv, scomplex *h, integer *ldh, scomplex *work)
{
  clasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}
inline void lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *h, integer *ldh, dcomplex *work)
{
  zlasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}

// --- computes a partial factorization of a float symmetric indefinite matrix using bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline void lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda, float *e, integer *ipiv, float *w, integer * ldw, integer *info)
{
  slasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline void lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda, double *e, integer *ipiv, double *w, integer * ldw, integer *info)
{
  dlasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline void lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, scomplex *e, integer *ipiv, scomplex *w, integer * ldw, integer *info)
{
  clasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline void lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, dcomplex *w, integer * ldw, integer *info)
{
  zlasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}

// --- computes a partial factorization of a float symmetric matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline void lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda, integer *ipiv, float *w, integer * ldw, integer *info)
{
  slasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline void lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda, integer *ipiv, double *w, integer * ldw, integer *info)
{
  dlasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline void lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer * ldw, integer *info)
{
  clasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline void lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer * ldw, integer *info)
{
  zlasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}

// --- solves a triangular banded system of equations ---
inline void latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, float *ab, integer *ldab, float *x, float *scale, float *cnorm, integer *info)
{
  slatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}
inline void latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, double *ab, integer *ldab, double *x, double *scale, double *cnorm, integer *info)
{
  dlatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}
inline void latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, scomplex *ab, integer *ldab, scomplex *x, float *scale, float *cnorm, integer *info)
{
  clatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}
inline void latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, dcomplex *ab, integer *ldab, dcomplex *x, double *scale, double *cnorm, integer *info)
{
  zlatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}

// --- uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate---
inline void latdf(integer *ijob, integer *n, float *z, integer * ldz, float *rhs, float *rdsum, float *rdscal, integer *ipiv, integer * jpiv)
{
  slatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}
inline void latdf(integer *ijob, integer *n, double *z, integer * ldz, double *rhs, double *rdsum, double *rdscal, integer *ipiv, integer * jpiv)
{
  dlatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}
inline void latdf(integer *ijob, integer *n, scomplex *z, integer * ldz, scomplex *rhs, float *rdsum, float *rdscal, integer *ipiv, integer * jpiv)
{
  clatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}
inline void latdf(integer *ijob, integer *n, dcomplex *z, integer * ldz, dcomplex *rhs, double *rdsum, double *rdscal, integer *ipiv, integer * jpiv)
{
  zlatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}

// --- solves a triangular system of equations with the matrix held in packed storage ---
inline void latps(char *uplo, char *trans, char *diag, char * normin, integer *n, float *ap, float *x, float *scale, float *cnorm, integer *info)
{
  slatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}
inline void latps(char *uplo, char *trans, char *diag, char * normin, integer *n, double *ap, double *x, double *scale, double *cnorm, integer *info)
{
  dlatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}
inline void latps(char *uplo, char *trans, char *diag, char * normin, integer *n, scomplex *ap, scomplex *x, float *scale, float *cnorm, integer *info)
{
  clatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}
inline void latps(char *uplo, char *trans, char *diag, char * normin, integer *n, dcomplex *ap, dcomplex *x, double *scale, double *cnorm, integer *info)
{
  zlatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}

// --- reduces the first nb rows and columns of a symmetric/Hermitian matrix A to float tridiagonal form by an orthogonal similarity transformation ---
inline void latrd(char *uplo, integer *n, integer *nb, float *a, integer *lda, float *e, float *tau, float *w, integer *ldw)
{
  slatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}
inline void latrd(char *uplo, integer *n, integer *nb, double *a, integer *lda, double *e, double *tau, double *w, integer *ldw)
{
  dlatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}
inline void latrd(char *uplo, integer *n, integer *nb, scomplex *a, integer *lda, float *e, scomplex *tau, scomplex *w, integer *ldw)
{
  clatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}
inline void latrd(char *uplo, integer *n, integer *nb, dcomplex *a, integer *lda, double *e, dcomplex *tau, dcomplex *w, integer *ldw)
{
  zlatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}

// --- solves a triangular system of equations with the scale factor set to prevent overflow ---
inline void latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, float *a, integer *lda, float *x, float *scale, float *cnorm, integer *info)
{
  slatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}
inline void latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, double *a, integer *lda, double *x, double *scale, double *cnorm, integer *info)
{
  dlatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}
inline void latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, scomplex *a, integer *lda, scomplex *x, float *scale, float *cnorm, integer *info)
{
  clatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}
inline void latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, dcomplex *a, integer *lda, dcomplex *x, double *scale, double *cnorm, integer *info)
{
  zlatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}

// --- factors an upper trapezoidal matrix by means of orthogonal transformations ---
inline void latrz(integer *m, integer *n, integer *l, float *a, integer *lda, float *tau, float *work)
{
  slatrz_(m, n, l, a, lda, tau, work);
}
inline void latrz(integer *m, integer *n, integer *l, double *a, integer *lda, double *tau, double *work)
{
  dlatrz_(m, n, l, a, lda, tau, work);
}
inline void latrz(integer *m, integer *n, integer *l, scomplex *a, integer *lda, scomplex *tau, scomplex *work)
{
  clatrz_(m, n, l, a, lda, tau, work);
}
inline void latrz(integer *m, integer *n, integer *l, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work)
{
  zlatrz_(m, n, l, a, lda, tau, work);
}

// --- computes a blocked Tall-Skinny QR factorization of a float M-by-N matrix A for M >= N ---
inline void latsqr(integer *m, integer *n, integer *mb, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *work, integer *lwork, integer *info)
{
  slatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void latsqr(integer *m, integer *n, integer *mb, integer *nb, double *a, integer *lda, double *t, integer *ldt, double *work, integer *lwork, integer *info)
{
  dlatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void latsqr(integer *m, integer *n, integer *mb, integer *nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  clatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void latsqr(integer *m, integer *n, integer *mb, integer *nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  zlatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}

inline float lantp(char* norm, char* uplo, char* diag, integer* n, float* ap, float* work)
{
  return slantp_(norm, uplo, diag, n, ap, work);
}
inline double lantp(char* norm, char* uplo, char* diag, integer* n, double* ap, double* work)
{
  return dlantp_(norm, uplo, diag, n, ap, work);
}
inline float lantp(char* norm, char* uplo, char* diag, integer* n, scomplex* ap, float* work)
{
  return clantp_(norm, uplo, diag, n, ap, work);
}
inline double lantp(char* norm, char* uplo, char* diag, integer* n, dcomplex* ap, double* work)
{
  return zlantp_(norm, uplo, diag, n, ap, work);
}

inline float lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, float* ab, integer* ldab, float* work)
{
  return slantb_(norm, uplo, diag, n, k, ab, ldab, work);
}
inline double lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, double* ab, integer* ldab, double* work)
{
  return dlantb_(norm, uplo, diag, n, k, ab, ldab, work);
}
inline float lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, scomplex* ab, integer* ldab, float* work)
{
  return clantb_(norm, uplo, diag, n, k, ab, ldab, work);
}
inline double lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, dcomplex* ab, integer* ldab, double* work)
{
  return zlantb_(norm, uplo, diag, n, k, ab, ldab, work);
}

inline void gelqt(integer* m, integer* n, integer* mb, float* a, integer* lda, float* t, integer* ldt, float* work, integer* info)
{
  sgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}
inline void gelqt(integer* m, integer* n, integer* mb, double* a, integer* lda, double* t, integer* ldt, double* work, integer* info)
{
  dgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}
inline void gelqt(integer* m, integer* n, integer* mb, scomplex* a, integer* lda, scomplex* t, integer* ldt, scomplex* work, integer* info)
{
  cgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}
inline void gelqt(integer* m, integer* n, integer* mb, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, dcomplex* work, integer* info)
{
  zgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}

inline void gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, float* v, integer* ldv, float* t, integer* ldt, float* c, integer* ldc, float* work, integer* info)
{
  sgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}
inline void gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, double* v, integer* ldv, double* t, integer* ldt, double* c, integer* ldc, double* work, integer* info)
{
  dgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}
inline void gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, scomplex* v, integer* ldv, scomplex* t, integer* ldt, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  cgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}
inline void gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, dcomplex* v, integer* ldv, dcomplex* t, integer* ldt, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  zgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}

inline void getc2(integer* n, float* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  sgetc2_(n, a, lda, ipiv, jpiv, info);
}
inline void getc2(integer* n, double* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  dgetc2_(n, a, lda, ipiv, jpiv, info);
}
inline void getc2(integer* n, scomplex* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  cgetc2_(n, a, lda, ipiv, jpiv, info);
}
inline void getc2(integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  zgetc2_(n, a, lda, ipiv, jpiv, info);
}

// --- LU factorization without pivoting ---
inline void getrfnp(integer* m, integer* n, float* a, integer* lda, integer* info)
{
  sgetrfnp_(m, n, a, lda, info);
}
inline void getrfnp(integer* m, integer* n, double* a, integer* lda, integer* info)
{
  dgetrfnp_(m, n, a, lda, info);
}
inline void getrfnp(integer* m, integer* n, scomplex* a, integer* lda, integer* info)
{
  cgetrfnp_(m, n, a, lda, info);
}
inline void getrfnp(integer* m, integer* n, dcomplex* a, integer* lda, integer* info)
{
  zgetrfnp_(m, n, a, lda, info);
}

inline void spffrt2(float  *ap, integer *n, integer * ncolm, float  *work, float  *work2)
{
  sspffrt2_(ap, n, ncolm, work, work2);
}
inline void spffrt2(double *ap, integer *n, integer * ncolm, double *work, double *work2)
{
  dspffrt2_(ap, n, ncolm, work, work2);
}
inline void spffrt2(scomplex *ap, integer *n, integer * ncolm, scomplex *work, scomplex *work2)
{
  cspffrt2_(ap, n, ncolm, work, work2);
}
inline void spffrt2(dcomplex *ap, integer *n, integer * ncolm, dcomplex *work, dcomplex *work2)
{
  zspffrt2_(ap, n, ncolm, work, work2);
}

inline void spffrtx(float *ap, integer *n, integer * ncolm, float *work, float *work2)
{
  sspffrtx_(ap, n, ncolm, work, work2);
}
inline void spffrtx(double *ap, integer *n, integer * ncolm, double *work, double *work2)
{
  dspffrtx_(ap, n, ncolm, work, work2);
}
inline void spffrtx(scomplex *ap, integer *n, integer * ncolm, scomplex *work, scomplex *work2)
{
  cspffrtx_(ap, n, ncolm, work, work2);
}
inline void spffrtx(dcomplex *ap, integer *n, integer * ncolm, dcomplex *work, dcomplex *work2)
{
  zspffrtx_(ap, n, ncolm, work, work2);
}

// --- LU factorization complete and incomplete without pivoting ---
inline void getrfnpi(integer *m, integer *n, integer *nfact, float *a, integer *lda, integer *info)
{
  sgetrfnpi_(m, n, nfact, a, lda, info);
}
inline void getrfnpi(integer *m, integer *n, integer *nfact, double *a, integer * lda, integer *info)
{
  dgetrfnpi_(m, n, nfact, a, lda, info);
}
inline void getrfnpi(integer *m, integer *n, integer *nfact, scomplex *a, integer *lda, integer *info)
{
  cgetrfnpi_(m, n, nfact, a, lda, info);
}
inline void getrfnpi(integer *m, integer *n, integer*nfact, dcomplex *a, integer *lda, integer *info)
{
  zgetrfnpi_(m, n, nfact, a, lda, info);
}

// NB2-sized column blocked QR-factorization
inline void getsqrhrt(integer *m, integer *n, integer *mb1, integer *nb1, integer *nb2, float *a, integer *lda, float *t, integer * ldt, float *work, integer *lwork, integer *info)
{
  sgetsqrhrt_(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info);
}
inline void getsqrhrt(integer *m, integer *n, integer *mb1, integer *nb1, integer *nb2, double *a, integer *lda, double * t, integer *ldt, double *work, integer *lwork, integer *info)
{
  dgetsqrhrt_(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info);
}
inline void getsqrhrt(integer *m, integer *n, integer *mb1, integer *nb1, integer *nb2, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  cgetsqrhrt_(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info);
}
inline void getsqrhrt(integer *m, integer *n, integer *mb1, integer *nb1, integer *nb2, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  zgetsqrhrt_(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info);
}

/* LAQZ0 computes the eigenvalues of a real matrix pair (H,T),
   where H is an upper Hessenberg matrix and T is upper triangular,
   using the double-shift QZ method.*/
inline void slaqz0(char *wants, char *wantq, char *wantz, integer * n, integer *ilo, integer *ihi, float *a, integer *lda, float *b, integer *ldb, float *alphar, float *alphai, float *beta, float *q, integer *ldq, float *z, integer *ldz, float *work, integer *lwork, integer *rec, integer *info)
{
  slaqz0_(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, rec, info);
}
inline void laqz0(char *wants, char *wantq, char *wantz, integer * n, integer *ilo, integer *ihi, double *a, integer *lda, double *b, integer *ldb, double *alphar, double *alphai, double *beta, double *q, integer *ldq, double *z, integer *ldz, double *work, integer *lwork, integer *rec, integer *info)
{
  dlaqz0_(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, rec, info);
}
inline void laqz0(char *wants, char *wantq, char *wantz, integer * n, integer *ilo, integer *ihi, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *alpha, scomplex *beta, scomplex *q, integer *ldq, scomplex *z, integer *ldz, scomplex *work, integer *lwork, float * rwork, integer *rec, integer *info)
{
  claqz0_(wants, wantq, wantz, n, ilo, ihi, a, lda,b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, rec, info);
}
inline void laqz0(char *wants, char *wantq, char *wantz, integer * n, integer *ilo, integer *ihi, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *alpha, dcomplex * beta, dcomplex *q, integer *ldq, dcomplex *z, integer * ldz, dcomplex *work, integer *lwork, double *rwork, integer * rec, integer *info)
{
  zlaqz0_(wants, wantq, wantz, n, ilo, ihi, a, lda,b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, rec, info);
}

/* Given a 3-by-3 matrix pencil (A,B), LAQZ1 sets v to a
   scalar multiple of the first column of the product.*/
inline void laqz1(float *a, integer *lda, float *b, integer *ldb, float *sr1, float *sr2, float *si, float *beta1, float *beta2, float *v)
{
  slaqz1_(a, lda, b, ldb, sr1, sr2, si, beta1, beta2, v);
}
inline void laqz1(double *a, integer *lda, double *b, integer *ldb, double *sr1, double *sr2, double *si, double *beta1, double *beta2, double *v)
{
  dlaqz1_(a, lda, b, ldb, sr1, sr2, si, beta1, beta2, v);
}
inline void laqz1(logical *ilq, logical *ilz, integer *k, integer * istartm, integer *istopm, integer *ihi, scomplex *a, integer *lda, scomplex *b, integer *ldb, integer *nq, integer *qstart, scomplex *q, integer *ldq, integer *nz, integer *zstart, scomplex *z, integer * ldz)
{
  claqz1_(ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz, zstart, z, ldz);
}
inline void laqz1(logical *ilq, logical *ilz, integer *k, integer * istartm, integer *istopm, integer *ihi, dcomplex *a, integer * lda, dcomplex *b, integer *ldb, integer *nq, integer *qstart, dcomplex *q, integer *ldq, integer *nz, integer *zstart, dcomplex *z, integer *ldz)
{
  zlaqz1_(ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz, zstart, z, ldz);
}

// LAQZ2 chases a 2x2 shift bulge in a matrix pencil down a single position.
inline void laqz2(logical *ilq, logical *ilz, integer *k, integer * istartm, integer *istopm, integer *ihi, float *a, integer *lda, float * b, integer *ldb, integer *nq, integer *qstart, float *q, integer *ldq, integer *nz, integer *zstart, float *z, integer *ldz)
{
  slaqz2_(ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz, zstart, z, ldz);
}
inline void laqz2(logical *ilq, logical *ilz, integer *k, integer * istartm, integer *istopm, integer *ihi, double *a, integer *lda, double *b, integer *ldb, integer *nq, integer *qstart, double *q, integer *ldq, integer *nz, integer *zstart, double *z, integer *ldz)
{
  dlaqz2_(ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz, zstart, z, ldz);
}
inline void laqz2(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nw, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *q, integer *ldq, scomplex *z, integer *ldz, integer *ns, integer *nd, scomplex *alpha, scomplex *beta, scomplex *qc, integer *ldqc, scomplex *zc, integer *ldzc, scomplex *work, integer *lwork, float *rwork, integer *rec, integer * info)
{
  claqz2_(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, ns, nd, alpha, beta, qc, ldqc, zc, ldzc, work, lwork, rwork, rec, info);
}
inline void laqz2(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nw, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *q, integer *ldq, dcomplex *z, integer *ldz, integer *ns, integer * nd, dcomplex *alpha, dcomplex *beta, dcomplex *qc, integer *ldqc, dcomplex *zc, integer *ldzc, dcomplex *work, integer *lwork, double *rwork, integer *rec, integer *info)
{
  zlaqz2_(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, ns, nd, alpha, beta, qc, ldqc, zc, ldzc, work, lwork, rwork, rec, info);
}

// LAQZ3 performs AED
inline void laqz3(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nw, float *a, integer *lda, float *b, integer *ldb, float *q, integer *ldq, float *z, integer *ldz, integer *ns, integer *nd, float *alphar, float *alphai, float *beta, float *qc, integer *ldqc, float *zc, integer *ldzc, float * work, integer *lwork, integer *rec, integer *info)
{
  slaqz3_(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, ns, nd, alphar, alphai, beta, qc, ldqc, zc, ldzc, work, lwork, rec, info);
}
inline void laqz3(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nw, double *a, integer *lda, double *b, integer *ldb, double *q, integer * ldq, double *z, integer *ldz, integer *ns, integer *nd, double *alphar, double *alphai, double *beta, double * qc, integer *ldqc, double *zc, integer *ldzc, double *work, integer *lwork, integer *rec, integer *info)
{
  dlaqz3_(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, ns, nd, alphar, alphai, beta, qc, ldqc, zc, ldzc, work, lwork, rec, info);
}
inline void laqz3(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nshifts, integer * nblock_desired, scomplex *alpha, scomplex *beta, scomplex *a, integer * lda, scomplex *b, integer *ldb, scomplex *q, integer *ldq, scomplex *z, integer *ldz, scomplex *qc, integer *ldqc, scomplex *zc, integer *ldzc, scomplex *work, integer *lwork, integer *info)
{
  claqz3_(ilschur, ilq, ilz, n, ilo, ihi, nshifts, nblock_desired, alpha, beta, a, lda, b, ldb, q, ldq, z, ldz, qc, ldqc, zc, ldzc, work, lwork, info);
}
inline void laqz3(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nshifts, integer * nblock_desired, dcomplex *alpha, dcomplex *beta, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *q, integer *ldq, dcomplex *z, integer *ldz, dcomplex *qc, integer *ldqc, dcomplex *zc, integer *ldzc, dcomplex *work, integer *lwork, integer *info)
{
  zlaqz3_(ilschur, ilq, ilz, n, ilo, ihi, nshifts, nblock_desired, alpha, beta, a, lda, b, ldb, q, ldq, z, ldz, qc, ldqc, zc, ldzc, work, lwork, info);
}

// LAQZ4 Executes a single multishift QZ sweep
inline void laqz4(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nshifts, integer * nblock_desired, float *sr, float *si, float *ss, float *a, integer *lda, float *b, integer *ldb, float *q, integer *ldq, float *z, integer * ldz, float *qc, integer *ldqc, float *zc, integer *ldzc, float *work, integer *lwork, integer *info)
{
  slaqz4_(ilschur, ilq, ilz, n, ilo, ihi, nshifts, nblock_desired, sr, si, ss, a, lda, b, ldb, q, ldq, z, ldz, qc, ldqc, zc, ldzc, work, lwork, info);
}
inline void laqz4(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nshifts, integer * nblock_desired, double *sr, double *si, double *ss, double *a, integer *lda, double *b, integer *ldb, double * q, integer *ldq, double *z, integer *ldz, double *qc, integer *ldqc, double *zc, integer *ldzc, double *work, integer *lwork, integer *info)
{
  dlaqz4_(ilschur, ilq, ilz, n, ilo, ihi, nshifts, nblock_desired, sr, si, ss, a, lda, b, ldb, q, ldq, z, ldz, qc, ldqc, zc, ldzc, work, lwork, info);
}

// LARFB_GETT applies a real Householder block reflector H from the left to a real (K+M)-by-N  "triangular-pentagonal" matrix
inline void larfb_gett(char *ident, integer *m, integer *n, integer *k, float *t, integer *ldt, float *a, integer *lda, float *b, integer *ldb, float *work, integer *ldwork)
{
  slarfb_gett_(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline void larfb_gett(char *ident, integer *m, integer *n, integer *k, double *t, integer *ldt, double *a, integer *lda, double *b, integer *ldb, double *work, integer *ldwork)
{
  dlarfb_gett_(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline void larfb_gett(char *ident, integer *m, integer *n, integer *k, scomplex *t, integer *ldt, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *work, integer *ldwork)
{
  clarfb_gett_(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline void larfb_gett(char *ident, integer *m, integer *n, integer *k, dcomplex *t, integer *ldt, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *work, integer * ldwork)
{
  zlarfb_gett_(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork);
}

// ORGTSQR_ROW generates an M-by-N real matrix Q_out with orthonormal columns from the output of LATSQR.
inline void gtsqr_row(integer *m, integer *n, integer *mb, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *work, integer *lwork, integer *info)
{
  sorgtsqr_row_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void gtsqr_row(integer *m, integer *n, integer *mb, integer *nb, double *a, integer *lda, double *t, integer *ldt, double *work, integer *lwork, integer *info)
{
  dorgtsqr_row_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void gtsqr_row(integer *m, integer *n, integer *mb, integer *nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  cungtsqr_row_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline void gtsqr_row(integer *m, integer *n, integer *mb, integer *nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  zungtsqr_row_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}

} // namespace libflame
#endif  //  #ifndef LIBFLAME_HH
