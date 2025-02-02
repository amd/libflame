#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define hetrf_rook_free() \
       free (a   ); \
       free (aref); \
       free (ipiv); \
       free (ipivref)


/* Begin hetrf_rook_scomplex_parameters  class definition */
class hetrf_rook_scomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      hetrf_rook_scomplex_parameters ( int matrix_layout_i,char uplo_i,
                                 lapack_int n_i,lapack_int lda_i);
      ~hetrf_rook_scomplex_parameters (); 
};  /* end of hetrf_rook_scomplex_parameters  class definition */


/* Constructor hetrf_rook_scomplex_parameters definition */
hetrf_rook_scomplex_parameters:: hetrf_rook_scomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                            lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       hetrf_rook_free();
       printf(" hetrf_rook_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

hetrf_rook_scomplex_parameters:: ~hetrf_rook_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrf_rook_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_rook_free();
}

TEST(hetrf_rook,chetrf_rook1) {

    /* LAPACKE CHETRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf_rook) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                                                      lapack_int *ipiv);

    Fptr_NL_LAPACKE_chetrf_rook CHETRF_ROOK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hetrf_rook_scomplex_parameters   chetrf_rook_obj(LAPACK_COL_MAJOR,'U',510,521);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHETRF_ROOK = (Fptr_NL_LAPACKE_chetrf_rook)dlsym(hModule,"LAPACKE_chetrf_rook");
    if (NULL == CHETRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chetrf_rook_obj.inforef = CHETRF_ROOK( chetrf_rook_obj.matrix_layout,chetrf_rook_obj.uplo,
                            chetrf_rook_obj.n,chetrf_rook_obj.aref,chetrf_rook_obj.lda,
                                                      chetrf_rook_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chetrf_rook_obj.info = LAPACKE_chetrf_rook(chetrf_rook_obj.matrix_layout,chetrf_rook_obj.uplo,
                                    chetrf_rook_obj.n,chetrf_rook_obj.a,chetrf_rook_obj.lda,
                                                             chetrf_rook_obj.ipiv );

    if( chetrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrf_rook is wrong\n",
                    chetrf_rook_obj.info );
    }
    if( chetrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf_rook is wrong\n",
        chetrf_rook_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chetrf_rook_obj.n,chetrf_rook_obj.ipiv,chetrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chetrf_rook_obj.n*chetrf_rook_obj.lda,chetrf_rook_obj.a,chetrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(hetrf_rook,chetrf_rook2) {

    /* LAPACKE CHETRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf_rook) ( int matrix_layout ,char uplo ,
                   lapack_int n ,lapack_complex_float *a ,lapack_int lda,
                                                      lapack_int *ipiv);

    Fptr_NL_LAPACKE_chetrf_rook CHETRF_ROOK;
    void *hModule,*dModule;
    float diff;
    int ipiv_diff;
    hetrf_rook_scomplex_parameters   chetrf_rook_obj(LAPACK_COL_MAJOR,'L',421,478);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    CHETRF_ROOK = (Fptr_NL_LAPACKE_chetrf_rook)dlsym(hModule,"LAPACKE_chetrf_rook");
    if (NULL == CHETRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chetrf_rook_obj.inforef = CHETRF_ROOK( chetrf_rook_obj.matrix_layout,chetrf_rook_obj.uplo,
                            chetrf_rook_obj.n,chetrf_rook_obj.aref,chetrf_rook_obj.lda,
                                                      chetrf_rook_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    chetrf_rook_obj.info = LAPACKE_chetrf_rook(chetrf_rook_obj.matrix_layout,chetrf_rook_obj.uplo,
                                    chetrf_rook_obj.n,chetrf_rook_obj.a,chetrf_rook_obj.lda,
                                                             chetrf_rook_obj.ipiv );

    if( chetrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrf_rook is wrong\n",
                    chetrf_rook_obj.info );
    }
    if( chetrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf_rook is wrong\n",
        chetrf_rook_obj.inforef );
    }

    ipiv_diff = computeDiff_i( chetrf_rook_obj.n,chetrf_rook_obj.ipiv,chetrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( chetrf_rook_obj.n*chetrf_rook_obj.lda,chetrf_rook_obj.a,chetrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin hetrf_rook_dcomplex_parameters  class definition */
class hetrf_rook_dcomplex_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // upper or lower triangular part of A is stored
      lapack_int n; // No of rows,Columns
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array ab contains the matrix A
      /* Output parameters */
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      hetrf_rook_dcomplex_parameters ( int matrix_layout_i,char uplo_i,
                                lapack_int n_i,lapack_int lda_i);
      ~hetrf_rook_dcomplex_parameters (); 
};  /* end of hetrf_rook_dcomplex_parameters  class definition */


/* Constructor hetrf_rook_dcomplex_parameters definition */
hetrf_rook_dcomplex_parameters:: hetrf_rook_dcomplex_parameters ( int matrix_layout_i,
                                       char uplo_i,lapack_int n_i,
                                       lapack_int lda_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (lda*n)); 
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if( (a==NULL) || (aref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       hetrf_rook_free();
       printf(" hetrf_rook_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,lda*n);

   } /* end of Constructor  */

hetrf_rook_dcomplex_parameters:: ~hetrf_rook_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrf_rook_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_rook_free();
}

TEST(hetrf_rook,zhetrf_rook1) {

    /* LAPACKE ZHETRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_rook) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhetrf_rook ZHETRF_ROOK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hetrf_rook_dcomplex_parameters   zhetrf_rook_obj(LAPACK_ROW_MAJOR,'U',100,121);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHETRF_ROOK = (Fptr_NL_LAPACKE_zhetrf_rook)dlsym(hModule,"LAPACKE_zhetrf_rook");
    if (NULL == ZHETRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhetrf_rook_obj.inforef = ZHETRF_ROOK( zhetrf_rook_obj.matrix_layout,zhetrf_rook_obj.uplo,
                            zhetrf_rook_obj.n,zhetrf_rook_obj.aref,zhetrf_rook_obj.lda,
                                                     zhetrf_rook_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhetrf_rook_obj.info = LAPACKE_zhetrf_rook(zhetrf_rook_obj.matrix_layout,zhetrf_rook_obj.uplo,
                                   zhetrf_rook_obj.n,zhetrf_rook_obj.a,zhetrf_rook_obj.lda,
                                                             zhetrf_rook_obj.ipiv);

    if( zhetrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrf_rook is wrong\n",
                    zhetrf_rook_obj.info );
    }
    if( zhetrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf_rook is wrong\n",
        zhetrf_rook_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhetrf_rook_obj.n,zhetrf_rook_obj.ipiv,zhetrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zhetrf_rook_obj.n*zhetrf_rook_obj.lda,zhetrf_rook_obj.a,zhetrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}


TEST(hetrf_rook,zhetrf_rook2) {

    /* LAPACKE ZHETRF_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_rook) ( int matrix_layout ,char uplo ,
                  lapack_int n ,lapack_complex_double *a ,lapack_int lda,
                        lapack_int *ipiv);

    Fptr_NL_LAPACKE_zhetrf_rook ZHETRF_ROOK;
    void *hModule,*dModule;
    double diff;
    int ipiv_diff;
    hetrf_rook_dcomplex_parameters   zhetrf_rook_obj(LAPACK_ROW_MAJOR,'L',1000,1210);

    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    ZHETRF_ROOK = (Fptr_NL_LAPACKE_zhetrf_rook)dlsym(hModule,"LAPACKE_zhetrf_rook");
    if (NULL == ZHETRF_ROOK)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhetrf_rook_obj.inforef = ZHETRF_ROOK( zhetrf_rook_obj.matrix_layout,zhetrf_rook_obj.uplo,
                            zhetrf_rook_obj.n,zhetrf_rook_obj.aref,zhetrf_rook_obj.lda,
                                                     zhetrf_rook_obj.ipivref);

        /* Compute libflame's Lapacke o/p  */
    zhetrf_rook_obj.info = LAPACKE_zhetrf_rook(zhetrf_rook_obj.matrix_layout,zhetrf_rook_obj.uplo,
                                   zhetrf_rook_obj.n,zhetrf_rook_obj.a,zhetrf_rook_obj.lda,
                                                             zhetrf_rook_obj.ipiv);

    if( zhetrf_rook_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrf_rook is wrong\n",
                    zhetrf_rook_obj.info );
    }
    if( zhetrf_rook_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf_rook is wrong\n",
        zhetrf_rook_obj.inforef );
    }

    ipiv_diff = computeDiff_i( zhetrf_rook_obj.n,zhetrf_rook_obj.ipiv,zhetrf_rook_obj.ipivref);
    if( ipiv_diff >0){
        printf("\n warning: ipivot computation in dhetrf_rook1 test case failed \n");
    }

    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zhetrf_rook_obj.n*zhetrf_rook_obj.lda,zhetrf_rook_obj.a,zhetrf_rook_obj.aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
