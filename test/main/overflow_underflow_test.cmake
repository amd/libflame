# Overflow and Underflow test cases for APIs.
set(OVERFLOW_UNDERFLOW_VALUES_TEST_CASES "gesvd sdcz A A 1 10 1 1 10 -1 1 --imatrix=O"
            "gesvd sdcz S O 10 10 10 10 10 -1 1 --imatrix=O"
            "gesvd sdcz N N 150 250 150 150 250 -1 1 --imatrix=O"
            "gesvd sdcz S S 400 200 400 400 200 -1 1 --imatrix=O"
            "gesvd sdcz A A 1 10 1 1 10 -1 1 --imatrix=U"
            "gesvd sdcz S S 10 10 10 10 10 -1 1 --imatrix=U"
            "gesvd sdcz N A 150 250 150 150 250 -1 1 --imatrix=U"
            "gesvd sdcz A S 400 200 400 400 400 -1 1 --imatrix=U"
            "geev sdcz N V 100 100 100 100 -1 1 --imatrix=O"
            "geev sdcz V N 500 500 500 500 -1 1 --imatrix=O"
            "geev sdcz V V 1000 1000 1000 1000 -1 1 --imatrix=O"
            "geev sdcz N V 100 100 100 100 -1 1 --imatrix=U"
            "geev sdcz V N 500 500 500 500 -1 1 --imatrix=U"
            "geev sdcz V V 1000 1000 1000 1000 -1 1 --imatrix=U"
            "gerq2 sdcz 100 30 200 1 --imatrix=U"
            "gerq2 sdcz 75 30 200 1 --imatrix=O"
            "gerq2 sdcz 200 100 200 1 --imatrix=O"
            "potrf sdcz L 10 10 1 --imatrix=U"
            "potrf sdcz U 10 10 1 --imatrix=U"
            "potrf sdcz L 10 10 1 --imatrix=O"
            "potrf sdcz U 10 10 1 --imatrix=O"
            "potrf sdcz L 155 500 1 --imatrix=U"
            "potrf sdcz U 155 500 1 --imatrix=U"
            "potrf sdcz L 155 500 1 --imatrix=O"
            "potrf sdcz U 155 500 1 --imatrix=O"
            "potrf sdcz L 553 600 1 --imatrix=U"
            "potrf sdcz U 553 600 1 --imatrix=U"
            "potrf sdcz L 553 600 1 --imatrix=O"
            "potrf sdcz U 553 600 1 --imatrix=O"
            "gelqf sdcz 5 10 10 -1 1 --imatrix=U"
            "gelqf sdcz 100 79 199 -1 1 --imatrix=U"
            "gelqf sdcz 100 50 100 -1 1 --imatrix=O"
            "gelqf sdcz 20 10 20 -1 1 --imatrix=O"
        )

foreach(ou_vals_test_cases IN LISTS OVERFLOW_UNDERFLOW_VALUES_TEST_CASES)
    string(REPLACE " " ";" COMMANDLINE_PARAMS ${ou_vals_test_cases})
    set(TEST_NAME OVERFLOW_UNDERFLOW_VALUES_TEST_CASE_${TEST_NUM} )
    add_test(${TEST_NAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME} ${COMMANDLINE_PARAMS})
    set_tests_properties(${TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL;No test was run, give valid arguments")
MATH(EXPR TEST_NUM "${TEST_NUM}+1")
endforeach()

