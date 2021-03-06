/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  POTRF computes the Cholesky factorization of a symmetric (hermitian)
  positive definite matrix A.

  INFO is INTEGER
  = 0: successful exit
  < 0: if INFO = -i, the i-th argument had an illegal value -  LAPACK_potrf_op_check
  > 0: if INFO = i, the leading minor of order i is not - FLA_Chol
  positive definite, and the factorization could not be
  completed.
*/

extern void DTL_Trace(
		    uint8 ui8LogLevel,
		    uint8 ui8LogType,
		    const int8 *pi8FileName,
		    const int8 *pi8FunctionName,
		    uint32 ui32LineNumber,
		    const int8 *pi8Message);

#define FLA_ENABLE_ALT_PATHS  0

#define LAPACK_potrf(prefix)                                    \
  int F77_ ## prefix ## potrf( char* uplo,                      \
                               integer*  n,                         \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               integer*  ldim_A,                    \
                               integer*  info )

#define LAPACK_potrf_body(prefix)                                                            \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);                                              \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                                     \
  FLA_Uplo     uplo_fla;                                                                     \
  FLA_Obj      A;                                                                            \
  FLA_Error    e_val = FLA_SUCCESS;                                                          \
  FLA_Error    init_result;                                                                  \
  FLA_Bool skip = FALSE;                                                                     \
                                                                                             \
  if(datatype == FLA_FLOAT && !FLA_ENABLE_ALT_PATHS ){ /*sp doesnt go to lf*/                \
  if( *n < FLA_POTRF_FLOAT_SMALL )                                                           \
  lapack_spotf2( uplo, n, buff_A, ldim_A,  info );                                           \
  else                                                                                       \
  lapack_spotrf( uplo, n, buff_A, ldim_A,  info );                                           \
  if( info != 0 ) skip = TRUE;                                                               \
  }                                                                                          \
  /*for sizes greater than 423 dp uses lf code*/                                             \
  else if ( datatype == FLA_DOUBLE && *n < FLA_POTRF_DOUBLE_MID && !FLA_ENABLE_ALT_PATHS ) { \
  if( *n < FLA_POTRF_DOUBLE_SMALL )                                                          \
  lapack_dpotf2( uplo, n, buff_A, ldim_A,  info );                                           \
  else                                                                                       \
  lapack_dpotrf( uplo, n, buff_A, ldim_A,  info );                                           \
  if( info != 0 ) skip = TRUE;                                                               \
  }                                                                                          \
  else {                                                                                     \
  FLA_Init_safe( &init_result );                                \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );        \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *n, *n, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  e_val = FLA_Chol( uplo_fla, A );                              \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
  }                                                             \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                \
  else if( skip != TRUE)      *info = 0;                        \
                                                                \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                  \
  return 0;

LAPACK_potrf(s)
{
    {
        LAPACK_RETURN_CHECK( spotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(s)
    }
}
LAPACK_potrf(d)
{
    {
        LAPACK_RETURN_CHECK( dpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(d)
    }
}
LAPACK_potrf(c)
{
    {
        LAPACK_RETURN_CHECK( cpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(c)
    }
}
LAPACK_potrf(z)
{
    {
        LAPACK_RETURN_CHECK( zpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(z)
    }
}

#define LAPACK_potf2(prefix)                                    \
  int F77_ ## prefix ## potf2( char* uplo,                      \
                               integer*  n,                         \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               integer*  ldim_A,                    \
                               integer*  info )

LAPACK_potf2(s)
{
    {
        LAPACK_RETURN_CHECK( spotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(s)
    }
}
LAPACK_potf2(d)
{
    {
        LAPACK_RETURN_CHECK( dpotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(d)
    }
}
LAPACK_potf2(c)
{
    {
        LAPACK_RETURN_CHECK( cpotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(c)
    }
}
LAPACK_potf2(z)
{
    {
        LAPACK_RETURN_CHECK( zpotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(z)
    }
}

#endif
