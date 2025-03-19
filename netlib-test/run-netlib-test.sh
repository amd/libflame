#!/bin/bash

echo
echo "Argument Values"
LAPACK_TEST_DIR=lapack-3.12
BLAS_LIB=libblis-mt.a
BLAS_LIB_PATH= 
LAPACK_LIB=libflame.a
LAPACK_LIB_PATH=
AOCLUTILS_LIB_PATH=
AOCLUTILS_LIB=libaoclutils.a
DTL_LIB=libaocldtl.a
DTL_LIB_PATH=
ILP64=0
DTL=0
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
         AOCLUTILS_LIB_PATH) AOCLUTILS_LIB_PATH=${DATA} ;; 
         DTL_LIB_PATH)       DTL_LIB_PATH=${DATA} ;; 
         DTL_LIB)            DTL_LIB=${DATA} ;;
         LAPACK_TEST_DIR)    LAPACK_TEST_DIR=${DATA} ;;     
         ILP64)              ILP64=${DATA} ;;   
         DTL)                DTL=${DATA} ;; 
		 GCOV)               GCOV=${DATA} ;;  
		 AOCL_LAPACK_SUMMARY)  AOCL_LAPACK_SUMMARY=${DATA} ;;
         *)   
   esac    
done

echo "BLAS_LIB_PATH = $BLAS_LIB_PATH"
echo "BLAS_LIB = $BLAS_LIB"
echo "LAPACK_LIB_PATH = $LAPACK_LIB_PATH"
echo "LAPACK_LIB = $LAPACK_LIB"
echo "AOCLUTILS_LIB_PATH = $AOCLUTILS_LIB_PATH"
echo "LAPACK_TEST_DIR = $LAPACK_TEST_DIR"
echo
echo
echo "**********************************"
echo

if [[ $BLAS_LIB_PATH  == "" || $LAPACK_LIB_PATH == "" || $AOCLUTILS_LIB_PATH == "" ]]
then
	echo "Error in calling script"
        echo "----------------------------------"
	echo
	echo "Usage :"
	echo
	echo "$ sh run-netlib-test.sh BLAS_LIB_PATH=<blas library path> LAPACK_LIB_PATH=<lapack library path> "
	echo "     AOCLUTILS_LIB_PATH=<aocl-utils library path> [BLAS_LIB=<blas library] [LAPACK_LIB=<lapack library>] "
	echo "     [ILP64=<0/1>] [LAPACK_TEST_DIR=<netlib lapack test directory name>]"
	echo "     [GCOV=<0/1>]"
	echo "[] indicates optional argument"
	echo 
	echo "Example: $ sh run-netlib-test.sh BLAS_LIB_PATH=\"/home/user/blis/install/lib\" LAPACK_LIB_PATH=\"/home/user/libflame/install/lib\" AOCLUTILS_LIB_PATH=\"/home/user/aoclutils/install/lib\" BLAS_LIB=\"libblis.a\" LAPACK_LIB=\"libflame.a\""
  	echo
  	echo "BLAS_LIB : blas library to use. Default=libblist-mt.a"
	echo "LAPACK_LIB : lapac library to use. Default=libflame.a"
  	echo "DTL_LIB : DTL library binary to be used. Default=libaocldtl.a" 
  	echo "DTL : Enable(1) or disable(0) DTL. Default=0"
	echo "ILP64 : LP64 or ILP64 mode. Default=0(Use LP64)"
	echo "BLAS_LIB_PATH : path of blas library chosen in BLAS_LIB"
	echo "LAPACK_LIB_PATH : path to lapack library chosen in LAPACK_LIB"
	echo "AOCLUTILS_LIB_PATH : path to aocl-utils library"
  	echo "DTL_LIB_PATH : path to DTL library chosen in DTL_LIB (if DTL is enabled)"
	echo "LAPACK_TEST_DIR : netlib lapack test directory name. Default=lapack-3.12.0"
	echo "GCOV : Enable(1) or disable(0) Code Coverage. Only Enable if Code Coverage is enabled on the library. Default=0"
	echo "AOCL_LAPACK_SUMMARY : run aocl-lapack netlib test suite summary script. Default=1"
	echo
	exit 1
fi

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

if [[ $DTL != "0" ]]
then
  cp $DTL_LIB_PATH/$DTL_LIB .
fi

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
TESTLAPACKLIB="$PWD/liblapack.a $AOCLUTILS_LIB_PATH/$AOCLUTILS_LIB"

if [[ $ILP64 =~ ("1"|"ON") ]]
then
	FORTRAN_FLAGS="gfortran -fopenmp -fdefault-integer-8"
fi

if [[ $DTL = "1" ]]
then
	TESTLAPACKLIB="$PWD/liblapack.a $PWD/libaocldtl.a -lpthread"
fi

if [[ $GCOV =~ ("1"|"ON") ]]
then
	# echo "======>COVERAGE=1"
	GCOV_FLAGS="-lgcov --coverage"
fi

OMP_NUM_THREADS=1 make FC="$FORTRAN_FLAGS" LDFLAGS+="-lstdc++ -lpthread -fopenmp $GCOV_FLAGS" LAPACKLIB="$TESTLAPACKLIB" -j

if [[ $AOCL_LAPACK_SUMMARY = "1" ]]
then
	cd ../
	python3 run-netlib-test-summary.py
fi
