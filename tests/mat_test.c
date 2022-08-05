#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "../src/matrix.h"
#include <stdio.h>
#include <time.h>

int main (void)
{
   clock_t start_t, end_t;
   double total_t;
   int i;

   start_t = clock();
   printf("Starting of the program, start_t = %ld\n", start_t);
    
   printf("Going to scan a big loop, start_t = %ld\n", start_t);
   for(i=0; i< 10000000; i++) {
   }
   end_t = clock();
   printf("End of the big loop, end_t = %ld\n", end_t);
   
   total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
   printf("Total time taken by CPU: %f\n", total_t  );
   printf("Exiting of the program...\n");

   return(0);
  // Py_Initialize(); // Need to call this so that Python.h functions won't segfault
  // CU_pSuite pSuite = NULL;

  // /* initialize the CUnit test registry */
  // if (CU_initialize_registry() != CUE_SUCCESS)
    // return CU_get_error();

  // /* add a suite to the registry */
  // pSuite = CU_add_suite("mat_test_suite", init_suite, clean_suite);
  // if (pSuite == NULL) {
    // CU_cleanup_registry();
    // return CU_get_error();
  // }

   // /* add the tests to the suite */
   // if ((CU_add_test(pSuite, "add_test", add_test) == NULL) ||
        // // (OPTIONAL) Uncomment the following lines if you have implemented sub_matrix and neg_matrix.
        // (CU_add_test(pSuite, "sub_test", sub_test) == NULL) ||
        // (CU_add_test(pSuite, "neg_test", neg_test) == NULL) ||
        
        // (CU_add_test(pSuite, "mul_square_test", mul_square_test) == NULL) ||
        // (CU_add_test(pSuite, "mul_non_square_test", mul_non_square_test) == NULL) ||
        // (CU_add_test(pSuite, "abs_test", abs_test) == NULL) ||
        // (CU_add_test(pSuite, "pow_test", pow_test) == NULL) ||
        // (CU_add_test(pSuite, "alloc_fail_test", alloc_fail_test) == NULL) ||
        // (CU_add_test(pSuite, "alloc_success_test", alloc_success_test) == NULL) ||
        // (CU_add_test(pSuite, "alloc_ref_fail_test", alloc_ref_fail_test) == NULL) ||
        // (CU_add_test(pSuite, "alloc_ref_success_test", alloc_ref_success_test) == NULL) ||
        // (CU_add_test(pSuite, "dealloc_null_test", dealloc_null_test) == NULL) ||
        // (CU_add_test(pSuite, "get_test", get_test) == NULL) ||
        // (CU_add_test(pSuite, "set_test", set_test) == NULL)
     // )
   // {
      // CU_cleanup_registry();
      // return CU_get_error();
   // }

  // // Run all tests using the basic interface
  // // CU_basic_set_mode(CU_BRM_NORMAL);
  // CU_basic_set_mode(CU_BRM_VERBOSE);
  // CU_basic_run_tests();
  // printf("\n");
  // CU_basic_show_failures(CU_get_failure_list());
  // printf("\n\n");

  // /* Clean up registry and return */
  // CU_cleanup_registry();
  // return CU_get_error();
}
