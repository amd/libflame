
/******************************************************************************
* Copyright (c) 2019 - present Advanced Micro Devices, Inc. All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*******************************************************************************/

/*! @file liblame_orgbr.cc
 *  libflame_orgbr.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void orgbr_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = fla_min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  char vect = 'P';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::orgbr( LAPACK_COL_MAJOR, &vect, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_form_V_external( aCIOObj, taupCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "orgbr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orgbr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] aCPPIOBuff ;
  delete[] tauqCPPOBuff ;
  delete[] taupCPPOBuff ;
  delete[] aCIOBuff ;
  delete[] tauqCOBuff ;
  delete[] taupCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &tauqCOObj );
  FLA_Obj_free_without_buffer( &taupCOObj );
}

void orgbr_testall_variants(){
  orgbr_test<float>();
  orgbr_test<double>();
}

int main(int argc, char *argv[])
{
  orgbr_testall_variants();
}