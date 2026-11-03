#!/bin/bash

#
# Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
#

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/netlib-test-common.sh"

echo
echo "Argument Values"
LAPACK_TEST_DIR=lapack-3.12
BLAS_LIB=
BLAS_LIB_PATH= 
LAPACK_LIB=
LAPACK_LIB_PATH=
ILP64=0
GCOV=0
AOCL_LAPACK_SUMMARY=1

for ARG in "$@"
do
   VAR=$(echo $ARG | cut -f1 -d=)
   DATA=$(echo $ARG | cut -f2 -d=)   

   case "$VAR" in
         BLAS_LIB)           BLAS_LIB=${DATA} ;;
         LAPACK_LIB)         LAPACK_LIB=${DATA} ;;     
         BLAS_LIB_PATH)      BLAS_LIB_PATH=${DATA} ;;
         LAPACK_LIB_PATH)    LAPACK_LIB_PATH=${DATA} ;; 
         LAPACK_TEST_DIR)    LAPACK_TEST_DIR=${DATA} ;;     
         ILP64)              ILP64=${DATA} ;;   
		 GCOV)               GCOV=${DATA} ;;  
		 AOCL_LAPACK_SUMMARY)  AOCL_LAPACK_SUMMARY=${DATA} ;;
         *)   
   esac    
done

if [[ $BLAS_LIB_PATH  == "" || $LAPACK_LIB_PATH == "" ]]
then
	echo "Error in calling script"
        echo "----------------------------------"
	echo
	echo "Usage :"
	echo
	echo "$ sh run-netlib-test.sh BLAS_LIB_PATH=<blas library path> LAPACK_LIB_PATH=<lapack library path> "
	echo "     [BLAS_LIB=<blas library] [LAPACK_LIB=<lapack library>] "
	echo "     [ILP64=<0/1>] [LAPACK_TEST_DIR=<netlib lapack test directory name>]"
	echo "     [GCOV=<0/1>]"
	echo "[] indicates optional argument"
	echo
	echo "Example: $ sh run-netlib-test.sh BLAS_LIB_PATH=\"/home/user/blis/install/lib\" LAPACK_LIB_PATH=\"/home/user/libflame/install/lib\" BLAS_LIB=\"libblis.a\" LAPACK_LIB=\"libflame.a\""
  	echo
  	echo "BLAS_LIB : BLAS library to use. Default=libblis-mt.so, then libblis-mt.a"
	echo "LAPACK_LIB : LAPACK library to use. Default=libflame.so, then libflame.a"
	echo "ILP64 : LP64 or ILP64 mode. Default=0(Use LP64)"
	echo "BLAS_LIB_PATH : path of blas library chosen in BLAS_LIB"
	echo "LAPACK_LIB_PATH : path to lapack library chosen in LAPACK_LIB"
	echo "LAPACK_TEST_DIR : netlib lapack test directory name. Default=lapack-3.12.0"
	echo "GCOV : Enable(1) or disable(0) Code Coverage. Only Enable if Code Coverage is enabled on the library. Default=0"
	echo "AOCL_LAPACK_SUMMARY : run aocl-lapack netlib test suite summary script. Default=1"
	echo
	exit 1
fi

configure_test_libraries || exit 1

if [ ! -d $LAPACK_TEST_DIR ]
then
	echo "Directory $LAPACK_TEST_DIR does not exist. Exiting..."
	exit 1
fi

rm -rf libflame_netlib

./create_new_testdir.sh $LAPACK_TEST_DIR libflame_netlib
cd libflame_netlib

cp $BLAS_LIB_PATH/$BLAS_LIB .
cp $LAPACK_LIB_PATH/$LAPACK_LIB .

if [[ $LAPACK_LIB != "liblapack.a" ]]
then
   ln -s $LAPACK_LIB liblapack.a
fi

if [[ $BLAS_LIB != "libblas.a" ]]
then
    ln -s $BLAS_LIB libblas.a
fi

ulimit -s unlimited

FORTRAN_FLAGS="gfortran -fopenmp"
C_FLAGS="-O2"
TESTLAPACKLIB="$PWD/liblapack.a"

if [[ $ILP64 =~ ("1"|"ON") ]]
then
	FORTRAN_FLAGS="gfortran -fopenmp -fdefault-integer-8"
	C_FLAGS="${C_FLAGS} -DINT_64BIT"
fi

if [[ $GCOV =~ ("1"|"ON") ]]
then
	# echo "======>COVERAGE=1"
	GCOV_FLAGS="-lgcov --coverage"
fi

OMP_NUM_THREADS=1 make FC="$FORTRAN_FLAGS" CFLAGS="$C_FLAGS" LDFLAGS+="-lpthread -fopenmp $GCOV_FLAGS" LAPACKLIB="$TESTLAPACKLIB" -j
MAKE_STATUS=$?

if [[ $AOCL_LAPACK_SUMMARY = "1" ]]
then
	cd ../
	python3 run-netlib-test-summary.py
fi

exit $MAKE_STATUS
