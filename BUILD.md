# STEPS TO COMPILE AOCL-LAPACK USING CMAKE
## 1. Generating AOCL-LAPACK library

Create a new build directory e.g. build1
        
    mkdir build1
    cd build1

Use the following command to configure project

#### GCC:

    With LP64
        cmake ../ -DCMAKE_INSTALL_PREFIX=<path>
    
    With ILP64
        cmake ../ -DENABLE_ILP64=ON -DCMAKE_INSTALL_PREFIX=<path>
	
    Note: Use -DCMAKE_C_COMPILER flag to set the compiler
            -DCMAKE_C_COMPILER=gcc OR
            export CC=gcc

#### AOCC:

     export CC=clang
     export CXX=clang++
     export FC=flang
     export FLIBS="-lflang"

     With LP64
         cmake ../ -DCMAKE_INSTALL_PREFIX=<path>
     
     With ILP64
         cmake ../ -DENABLE_ILP64=ON -DCMAKE_INSTALL_PREFIX=<path>

    Note:
    Complex return ABI is now selected automatically based on the compiler (see "Complex return ABI" below).

#### AOC:

     export CC=aoc
     export CXX=aoc++
     export FC=aof

     With LP64
         cmake ../ -DCMAKE_INSTALL_PREFIX=<path>
     
     With ILP64
         cmake ../ -DENABLE_ILP64=ON -DCMAKE_INSTALL_PREFIX=<path>

#### Complex return ABI:

ENABLE_F2C_DOTC and ENABLE_VOID_RETURN_COMPLEX_FUNCTION control whether complex valued functions
(BLAS cdotc/cdotu/zdotc/zdotu and libFLAME routines such as cladiv/zladiv) return their result
through a hidden pointer argument instead of a complex return type. Their default is derived from
the configured C compiler:

| Platform | Compiler           | Default |
|----------|--------------------|---------|
| Windows  | All                | ON      |
| Linux    | clang/flang (AOCC) | ON      |
| Linux    | aoc/aof (AOCC)     | OFF     |
| Linux    | gcc                | OFF     |

Note:
ENABLE_AMD_FLAGS is ON by default, so the AMD optimized code paths are enabled without 
passing any extra option. Pass -DENABLE_AMD_FLAGS=OFF to build the reference (non AMD optimized) library.

Shared library is turned on by default. To generate Static library provide additional option

    -DBUILD_SHARED_LIBS=OFF

Compile library using following command. This will generate libflame.a/libflame.so library in the lib directory
        
    cmake --build . -j OR make -j

Install the library

    make install

Linking with AOCL-BLAS
------------------------------------
AOCL-LAPACK can be linked with any Netlib BLAS compliant library when compiled with standard cmake options above. However, AOCL-LAPACK provides an option to explicitly to link with AOCL-BLAS library at compile time. This option helps achieve better performance for certain APIs on AMD "Zen" CPUs by invoking lower level AOCL-BLAS APIs directly. This is controlled by ENABLE_AOCL_BLAS, which is ON by default, so a path to AOCL-BLAS must be available at configure time. Pass `-DENABLE_AOCL_BLAS=OFF` to build without the tight coupling

The path of AOCL-BLAS library can be provided in one of the following methods
1. Set "AOCL_ROOT" environment variable to the root path where AOCL-BLAS library `($AOCL_ROOT/lib)` and header files `("$AOCL_ROOT"/include)` are located. AOCL_ROOT must ideally be path where all AOCL libraries are installed including AOCL-Utils and AOCL-BLAS.
`export AOCL_ROOT=<path to AOCL-libraries>`

2. Specify root path of AOCL-BLAS library through cmake option "AOCL_ROOT"
`cmake -DAOCL_ROOT=<path to AOCL-BLAS> ...`

The path specified in AOCL_ROOT must have "include" directory and a "lib" directory that contains the necesaary header files and AOCL-BLAS binary respectively.

AOCL-Utils CPUID header
-----------------------
AOCL-LAPACK uses the header-only AOCL-Utils CPUID implementation. It is
inlined into libflame at compile time, so no AOCL-Utils library path or
binary library is required. The header can be supplied with
`-DLIBAOCLUTILS_INCLUDE_PATH=<path to the AOCL-Utils SDK/Include directory>`.


## 2. Building main Test and AOCL_FLA_PROGRESS Test Suite
In order to build tests an additional flag, BUILD_TEST, must be set ON

    -DBUILD_TEST=ON 
    -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH="<path to blas library>" 
    -DEXT_BLAS_LIBNAME="<blas_lib_name>"
    -DBLAS_HEADER_PATH="<path to AOCL-BLAS header file blis.h>"

This will enable aocl progress feature tests and main test suite. It will generate test_libFLAME_aocl , test_lapack.x executables in the respective directories.

#### Note:
1. Building tests require an external BLAS library. Refer to the Readme in the respective test suite directory for more details.
2. EXT_BLAS_LIBNAME flag can be set with space separated binaries for BLAS libraries which may have more than one binary file. But all of them
   must be under the same directory as set in flag CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH.

Recomended to use AOCL-BLAS sharedlib with AOCL-LAPACK sharedlib

## 3 Building Legacy test and Netlib test
#### 1. To build Legacy test suite use
    -DBUILD_LEGACY_TEST=ON 
    -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH="<path to blas library>" 
    -DEXT_BLAS_LIBNAME=blas_lib_name
    -DBLAS_HEADER_PATH="<path to AOCL-BLAS header file blis.h>"

Note: On Windows, to build and run legacy test suite, a separate macro flag is enabled during AOCL-LAPACK library build because of certain constraints in legacy test suite.
    
#### 2. To Build Netlib-test add `-DBUILD_NETLIB_TEST=ON` along with cmake commands.
        
Note: Windows requires running create_new_testdir.bat script before running netlib test

## 4. Building documentation
    # Requirements:
        1. python >= 3.9
        2. doxygen >= 1.9.6
        3. sphinx >= 7.3.7
            python -m pip install sphinx==7.3.7
        4. rocm-docs-core
            python -m pip install rocm-docs-core
        5. breathe >= 4.30.0
            python -m pip install breathe==4.30.0

    # To build the documentation in /libflame/docs/libflame/sphinx/html directory use
     -DBUILD_DOC=ON
    e.g.
        cmake ../ -DBUILD_DOC=ON -DCMAKE_INSTALL_PREFIX=<path>

## 5. ENABLE TRACE and LOGS
    
User may also enable trace and logs by passing `-DENABLE_AOCL_DTL=[OPTION]`
along with setting the value of Macros `LF_AOCL_DTL_TRACE_ENABLE` and `LF_AOCL_DTL_LOG_ENABLE` to 1 in file `libflame/src/aocl_dtl/aocldtlcf.h`

e.g.
    
    cmake ../ -DENABLE_ILP64=OFF -DBUILD_TEST=ON -DENABLE_AOCL_DTL=[DTL_OPTION] -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH="<path to blas library>" -DEXT_BLAS_LIBNAME=<BLAS_lib_name> -DCMAKE_INSTALL_PREFIX=<path> -DBLAS_HEADER_PATH="<path to AOCL-BLAS header file blis.h>"

#### DTL_OPTIONS:

    1. "ALL" to ENABLE TRACE and LOG
	2. "TRACE" to ENABLE TRACE
	3. "LOG" to ENABLE LOGS
	4. "OFF" to Disable trace and log, if -DENABLE_AOCL_DTL is not passed with the cmake command, DTL is turned off

## 6. Using an external Lapack library to run tests

In order to run tests on an external lapack library an additional option
`-DEXT_LAPACK_LIBRARY_PATH="path/to/external/lapack/library"` and `-DEXT_LAPACK_LIBNAME="NAME_OF_THE_LAPACK_LIB"` can be passed, if the above options are left blank AOCL-LAPACK library will be used

## 7. Linking with an external openmp library
In order to link with an external openmp library user can pass `-DEXT_OPENMP_PATH=<openmp lib path> -DEXT_OPENMP_LIB=<openmp lib name>`
    
#### Note:   
1. In order to use openmp from the system -DEXT_OPENMP_PATH is to be left blank
2. To link Intel OpenMP library,libiomp5.so, set following flag addtionally
    - gcc: `-DCMAKE_C_FLAG="-fopenmp"`
    - aocc: `-DCMAKE_C_FLAG="-fopenmp=libiomp5"`


## 8. Using ctest
Ctest is enabled when `-DBUILD_TEST=ON` OR `-DBUILD_LEGACY_TEST=ON` OR `-DBUILD_NETLIB_TEST=ON`
To run ALL ctests together following command can be given.
    
    ctest --test-dir [BUILD_DIR]

To run a specific ctest following command can be given.
    
    ctest -R [TEST_NAME]

**Note**: Test names can be listed by
        
    ctest -N

To run build from any location
    
    ctest --test-dir [BUILD_DIR]

Additionally `--verbose` can be added to print the output from the executable.
#### Example:
Following command can be used to run tests with regular expression neg_test
            
    ctest --test-dir <build_dir> -R neg_test --verbose
    
On windows additional "-C Release"  is needed to run the test
    
    ctest --test-dir <build_dir> -R neg_test -C Release --verbose

To list all the tests `ctest --test-dir [BUILD_DIR] -N` can be given

## 9. ENABLE GCOV
In order to enable code coverage `-DENABLE_GCOV` can be passed during configuration.
After running the executable in the root directory run

    bash generate_code_coverage_html.sh

It will give you a prompt to view the code coverage of that particular application.

## 10. pkg-config support
When building the package using cmake, pkg-config metadata file for the library is generated in
the location `{CMAKE_INSTALL_PREFIX}/lib/pkgconfig/`. 
- The pkg-config metadata file should be placed at a location recognized by pkg-config tool.
- aocl-blas should be installed with its pkg-config metadata file. AOCL-Utils
  is consumed through its header-only CPUID interface.

## 11. Auto ISA detection
LF_ISA_CONFIG can be set to none, avx, avx2, avx2-strict, avx512, avx512-strict, or auto. When `-DLF_ISA_CONFIG=auto`, CMake runs a configure-time CPUID helper compiled with the configured C compiler
(CMAKE_C_COMPILER or CC). The helper checks AVX512F and AVX2 feature bits and selects avx512, avx2,
or none.

## Building using CMake presets

AOCL-LAPACK provides with build, testing and test workflow cmake presets for most common configurations. You can list
all available configurations with the following command:


| Preset type        |Command                                       |
|--------------------|----------------------------------------------|
| Configuration      | `cmake --list-presets`                       |
| Build              | `cmake --build --list-presets`               |
| CTest**            | `ctest --list-presets`                       |
| Workflow**         | `cmake --workflow --list-presets`            |

> ** CTest and Workflow presets are not available for Windows platform.

### Configuration presets

These presets are used to configure the project. AOCL-LAPACK requires the
AOCL-Utils CPUID header and requires a BLAS library when building tests.
The CPUID header is header-only and does not require an AOCL-Utils library.

> **Note for Windows:**
 On Windows, provide `LIBAOCLUTILS_INCLUDE_PATH` and `EXT_BLAS_LIBNAME` when
 configuring without the default header source and BLAS configuration.


Presets names use the following convention: `<os>-<buildsystem>-<compiler>-<st/mt>-<lp/ilp>-<static/shared>-<isa-mode>-<other optional commands>`

| Option                  | Allowed Values                                   | Notes                                                                                       |
|-------------------------|--------------------------------------------------|---------------------------------------------------------------------------------------------|
| os                      | linux, win (Windows)                             |                                                                                             |
| buildsystem             | - **Linux**: make <br> - **Windows**: ninja      | When using ninja on Windows, the build system will automatically use LLVM toolchain. **compiler** option should be omitted in this case |
| compiler                | - **Linux**: gcc, aocc <br> - **Windows**: msvc  | When using msvc, Microsoft Visual Studio toolchain is used as buildsystem. **buildsystem** option should be omitted in this case |
| st/mt                   | st, mt                                           | Single-threaded or multi-threaded build                                                     |
| lp/ilp                  | lp, ilp                                          | LP64 or ILP64 build                                                                         |
| static/shared           | static, shared                                   | Static or shared library                                                                    |
| isa-mode                | autoisa, avx, avx2, avx512                       | ISA mode for the build                                                                      |
| other optional commands | test, aoclblas                                   | - **test**: Builds libflame test binary <br> - **aoclblas**: Enables tight coupling with aoclblas |


Examples:

Building AOCL-LAPACK without tests:

- Using the default AOCL-Utils FetchContent header source

        cmake --preset <preset-name>

- Without the default AOCL-Utils header source

        cmake --preset <preset-name> -DLIBAOCLUTILS_INCLUDE_PATH=<path/to/libaoclutils/header/files>


Building AOCL-LAPACK with tests:

- Using the default AOCL-Utils header source with BLAS configuration

        cmake --preset <preset-name ending with test>
    
- Without the default AOCL-Utils header source and BLAS configuration
            
        cmake --preset <preset-name ending with test> -DLIBAOCLUTILS_INCLUDE_PATH=<path/to/libaoclutils/header/files> -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH=<path/to/blas/lib> -DEXT_BLAS_LIBNAME=<blas_lib_name> -DBLAS_HEADER_PATH=<path/to/blis.h>

> Note that when building libflame integrated with AOCL-BLAS library and if the blis pkg-config file is not available, AOCL-ROOT variable must be set as mentioned in the section "Linking with AOCL-BLAS".


### Build presets
AOCL-LAPACK provides two types of build presets for each configuration preset.

1. `build` preset: This preset is named same as the configuration preset. It is used to build the library.
1. `install` preset: This preset is named `install-<configuration-preset>`. It will build and install the library. By default, the library is installed in the `install-<suffix>` where suffix is determined by the configuration preset. If you want to install the library in a different location, you can provide the path using the `CMAKE_INSTALL_PREFIX` variable. while configuring the project.

```
cmake --build --preset <build-preset>
```

### Test presets
AOCL-LAPACK provides test presets for each configuration preset. These presets are used to run tests on the library. Once the library is built, you can run the tests using the following command:

    ctest --preset <test-preset>

### Workflow presets
Workflow presets are used to configure, build, install and run tests in a single command. You can use the following command to run the workflow:

    cmake --workflow --preset <workflow-preset>

As of now cmake workflow presets do not accept any additional arguments, so you cannot manually set the path to the AOCL-Utils header or BLAS library. Either use the pkg-config file for BLAS or as a workaround you can do the following:

1. Configure the project using the configuration preset with the required paths.
    
        cmake --preset <preset with same name as workflow preset> ... <additional arguments>

2. Run the workflow preset.
    
        cmake --workflow --preset <workflow-preset>

Note:
    autoconf based build system using configure and Makefile is no longer supported.
