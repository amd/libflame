#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sysv_aa_free() \
       if (ipiv != NULL) free (ipiv); \
       if (bref != NULL) free (bref); \
       if (b != NULL)    free (b   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (ipivref != NULL)free (ipivref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin sysv_aa_double_parameters  class definition */
class sysv_aa_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_aa_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_aa_double_parameters ();
};  /* end of sysv_aa_double_parameters  class definition */


/* Constructor sysv_aa_double_parameters definition */
sysv_aa_double_parameters:: sysv_aa_double_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_aa Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif


    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_aa_double_parameters object: malloc error.";
       sysv_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_aa_double_parameters:: ~sysv_aa_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_aa_double_parameters object: destructor invoked. \n");
#endif
   sysv_aa_free();
}


//  Test fixture class definition
class dsysv_aa_test  : public  ::testing::Test {
public:
   sysv_aa_double_parameters  *dsysv_aa_obj;
   void SetUp();  
   void TearDown () { delete dsysv_aa_obj; }
};


void dsysv_aa_test::SetUp(){

    /* LAPACKE DSYSV_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_dsysv_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsysv_aa DSYSV_AA;

    dsysv_aa_obj = new  sysv_aa_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsysv_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsysv_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsysv_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsysv_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYSV_AA = (Fptr_NL_LAPACKE_dsysv_aa)dlsym(dsysv_aa_obj->hModule, "LAPACKE_dsysv_aa");
    ASSERT_TRUE(DSYSV_AA != NULL) << "failed to get the Netlib LAPACKE_dsysv_aa symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    dsysv_aa_obj->inforef = DSYSV_AA( dsysv_aa_obj->matrix_layout,
                                  dsysv_aa_obj->uplo, dsysv_aa_obj->n,
                                  dsysv_aa_obj->nrhs,
                                  ( double *)dsysv_aa_obj->aref,
                              dsysv_aa_obj->lda, dsysv_aa_obj->ipivref,
                                dsysv_aa_obj->bref, dsysv_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsysv_aa_obj->info = LAPACKE_dsysv_aa( dsysv_aa_obj->matrix_layout,
                dsysv_aa_obj->uplo, dsysv_aa_obj->n, dsysv_aa_obj->nrhs,
                                  ( double *)dsysv_aa_obj->a,
                               dsysv_aa_obj->lda, dsysv_aa_obj->ipiv,
                                 dsysv_aa_obj->b, dsysv_aa_obj->ldb );


    if( dsysv_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsysv_aa is wrong\n",
                    dsysv_aa_obj->info );
    }
    if( dsysv_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsysv_aa is wrong\n",
        dsysv_aa_obj->inforef );
    }
}

TEST_F(dsysv_aa_test, dsysv_aa1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_aa_obj->b_bufsize,
                           dsysv_aa_obj->b, dsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_aa_test, dsysv_aa2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_aa_obj->b_bufsize,
                           dsysv_aa_obj->b, dsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_aa_test, dsysv_aa3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_aa_obj->b_bufsize,
                           dsysv_aa_obj->b, dsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_aa_test, dsysv_aa4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_aa_obj->b_bufsize,
                           dsysv_aa_obj->b, dsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sysv_aa_float_parameters  class definition */
class sysv_aa_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_aa_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sysv_aa_float_parameters ();
};  /* end of sysv_aa_float_parameters  class definition */


/* Constructor sysv_aa_float_parameters definition */
sysv_aa_float_parameters:: sysv_aa_float_parameters ( int matrix_layout_i,
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_aa float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_aa_double_parameters object: malloc error.";
       sysv_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sysv_aa_float_parameters:: ~sysv_aa_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_aa_float_parameters object: destructor invoked. \n");
#endif
   sysv_aa_free();
}


//  Test fixture class definition
class ssysv_aa_test  : public  ::testing::Test {
public:
   sysv_aa_float_parameters  *ssysv_aa_obj;
   void SetUp();  
   void TearDown () { delete ssysv_aa_obj; }
};


void ssysv_aa_test::SetUp(){

    /* LAPACKE SSYSV_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_ssysv_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  float * a,
                                  lapack_int lda,  lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssysv_aa SSYSV_AA;

    ssysv_aa_obj = new  sysv_aa_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssysv_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssysv_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssysv_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssysv_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYSV_AA = (Fptr_NL_LAPACKE_ssysv_aa)dlsym(ssysv_aa_obj->hModule, "LAPACKE_ssysv_aa");
    ASSERT_TRUE(SSYSV_AA != NULL) << "failed to get the Netlib LAPACKE_ssysv_aa symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    ssysv_aa_obj->inforef = SSYSV_AA( ssysv_aa_obj->matrix_layout,
                                  ssysv_aa_obj->uplo, ssysv_aa_obj->n,
                                  ssysv_aa_obj->nrhs,
                                  ( float *)ssysv_aa_obj->aref,
                                  ssysv_aa_obj->lda, ssysv_aa_obj->ipivref,
                                  ssysv_aa_obj->bref, ssysv_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssysv_aa_obj->info = LAPACKE_ssysv_aa( ssysv_aa_obj->matrix_layout,
                ssysv_aa_obj->uplo, ssysv_aa_obj->n, ssysv_aa_obj->nrhs,
                                  ( float *)ssysv_aa_obj->a,
                               ssysv_aa_obj->lda, ssysv_aa_obj->ipiv,
                                 ssysv_aa_obj->b, ssysv_aa_obj->ldb );
    if( ssysv_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssysv_aa is wrong\n",
                    ssysv_aa_obj->info );
    }
    if( ssysv_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssysv_aa is wrong\n",
        ssysv_aa_obj->inforef );
    }
}

TEST_F(ssysv_aa_test, ssysv_aa1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_aa_obj->b_bufsize,
                           ssysv_aa_obj->b, ssysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_aa_test, ssysv_aa2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_aa_obj->b_bufsize,
                           ssysv_aa_obj->b, ssysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_aa_test, ssysv_aa3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_aa_obj->b_bufsize,
                           ssysv_aa_obj->b, ssysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_aa_test, ssysv_aa4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_aa_obj->b_bufsize,
                           ssysv_aa_obj->b, ssysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_aa_scomplex_parameters  class definition */
class sysv_aa_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_aa_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_aa_scomplex_parameters ();
};  /* end of sysv_aa_scomplex_parameters  class definition */


/* Constructor sysv_aa_scomplex_parameters definition */
sysv_aa_scomplex_parameters:: sysv_aa_scomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_aa scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_aa_scomplex_parameters object: malloc error.";
       sysv_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_aa_scomplex_parameters:: ~sysv_aa_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_aa_scomplex_parameters object: destructor invoked. \n");
#endif
   sysv_aa_free();
}


//  Test fixture class definition
class csysv_aa_test  : public  ::testing::Test {
public:
   sysv_aa_scomplex_parameters  *csysv_aa_obj;
   void SetUp();  
   void TearDown () { delete csysv_aa_obj; }
};


void csysv_aa_test::SetUp(){

    /* LAPACKE CSYSV_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_csysv_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_float * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csysv_aa CSYSV_AA;

    csysv_aa_obj = new  sysv_aa_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csysv_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csysv_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csysv_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csysv_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYSV_AA = (Fptr_NL_LAPACKE_csysv_aa)dlsym(csysv_aa_obj->hModule, "LAPACKE_csysv_aa");
    ASSERT_TRUE(CSYSV_AA != NULL) << "failed to get the Netlib LAPACKE_csysv_aa symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    csysv_aa_obj->inforef = CSYSV_AA( csysv_aa_obj->matrix_layout,
                                  csysv_aa_obj->uplo, csysv_aa_obj->n,
                                  csysv_aa_obj->nrhs,
                                  ( lapack_complex_float *)csysv_aa_obj->aref,
                                  csysv_aa_obj->lda, csysv_aa_obj->ipivref,
                                  csysv_aa_obj->bref, csysv_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csysv_aa_obj->info = LAPACKE_csysv_aa( csysv_aa_obj->matrix_layout,
                csysv_aa_obj->uplo, csysv_aa_obj->n, csysv_aa_obj->nrhs,
                                  ( lapack_complex_float *)csysv_aa_obj->a,
                               csysv_aa_obj->lda, csysv_aa_obj->ipiv,
                                 csysv_aa_obj->b, csysv_aa_obj->ldb );


    if( csysv_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csysv_aa is wrong\n",
                    csysv_aa_obj->info );
    }
    if( csysv_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csysv_aa is wrong\n",
        csysv_aa_obj->inforef );
    }
}

TEST_F(csysv_aa_test, csysv_aa1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_aa_obj->b_bufsize,
                           csysv_aa_obj->b, csysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_aa_test, csysv_aa2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_aa_obj->b_bufsize,
                           csysv_aa_obj->b, csysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_aa_test, csysv_aa3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_aa_obj->b_bufsize,
                           csysv_aa_obj->b, csysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_aa_test, csysv_aa4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_aa_obj->b_bufsize,
                           csysv_aa_obj->b, csysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_aa_dcomplex_parameters  class definition */
class sysv_aa_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_aa_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_aa_dcomplex_parameters ();
};  /* end of sysv_aa_dcomplex_parameters  class definition */


/* Constructor sysv_aa_dcomplex_parameters definition */
sysv_aa_dcomplex_parameters:: sysv_aa_dcomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_aa DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_aa_dcomplex_parameters object: malloc error.";
       sysv_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_aa_dcomplex_parameters:: ~sysv_aa_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_aa_dcomplex_parameters object: destructor invoked. \n");
#endif
   sysv_aa_free();
}


//  Test fixture class definition
class zsysv_aa_test  : public  ::testing::Test {
public:
   sysv_aa_dcomplex_parameters  *zsysv_aa_obj;
   void SetUp();  
   void TearDown () { delete zsysv_aa_obj; }
};


void zsysv_aa_test::SetUp(){

    /* LAPACKE ZSYSV_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_zsysv_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsysv_aa ZSYSV_AA;

    zsysv_aa_obj = new  sysv_aa_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsysv_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsysv_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsysv_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsysv_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYSV_AA = (Fptr_NL_LAPACKE_zsysv_aa)dlsym(zsysv_aa_obj->hModule, "LAPACKE_zsysv_aa");
    ASSERT_TRUE(ZSYSV_AA != NULL) << "failed to get the Netlib LAPACKE_zsysv_aa symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zsysv_aa_obj->inforef = ZSYSV_AA( zsysv_aa_obj->matrix_layout,
                                  zsysv_aa_obj->uplo, zsysv_aa_obj->n,
                                  zsysv_aa_obj->nrhs,
                                  ( lapack_complex_double *)zsysv_aa_obj->aref,
                                  zsysv_aa_obj->lda, zsysv_aa_obj->ipivref,
                                  zsysv_aa_obj->bref, zsysv_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsysv_aa_obj->info = LAPACKE_zsysv_aa( zsysv_aa_obj->matrix_layout,
                zsysv_aa_obj->uplo, zsysv_aa_obj->n, zsysv_aa_obj->nrhs,
                                  ( lapack_complex_double *)zsysv_aa_obj->a,
                               zsysv_aa_obj->lda, zsysv_aa_obj->ipiv,
                                 zsysv_aa_obj->b, zsysv_aa_obj->ldb );


    if( zsysv_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsysv_aa is wrong\n",
                    zsysv_aa_obj->info );
    }
    if( zsysv_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsysv_aa is wrong\n",
        zsysv_aa_obj->inforef );
    }
}

TEST_F(zsysv_aa_test, zsysv_aa1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_aa_obj->b_bufsize,
                           zsysv_aa_obj->b, zsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_aa_test, zsysv_aa2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_aa_obj->b_bufsize,
                           zsysv_aa_obj->b, zsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_aa_test, zsysv_aa3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_aa_obj->b_bufsize,
                           zsysv_aa_obj->b, zsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_aa_test, zsysv_aa4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_aa_obj->b_bufsize,
                           zsysv_aa_obj->b, zsysv_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
