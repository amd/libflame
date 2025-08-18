###############################################################################
# Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
set(GROUND_TRUTH_TEST_CASES 
            "labrd scdz 75 110 50 150 150 150 10 --seed=100 --BRT=G"
            "labrd scdz 10 11 5 20 20 20 10 --seed=100 --BRT=G"
            "labrd scdz 80 50 40 100 100 100 10 --seed=100 --BRT=G"
            "labrd scdz 100 110 50 200 200 200 10 --seed=100 --BRT=G"
            "labrd scdz 110 100 75 150 150 150 10 --seed=100 --BRT=G"
            "labrd scdz 10 8 5 20 20 20 10 --seed=100 --BRT=G"
            "labrd scdz 100 180 80 200 200 200 10 --seed=100 --BRT=G"
            "labrd scdz 190 170 150 200 200 200 10 --seed=100 --BRT=G"
            "labrd scdz 200 190 150 200 200 200 10 --seed=100 --BRT=G"
            "labrd scdz 80 70 50 100 100 100 10 --seed=100 --BRT=G"
            "labrd scdz 190 170 150 200 200 200 10 --imatrix=U --seed=100 --BRT=G"
            "labrd scdz 200 190 150 200 200 200 10 --imatrix=O --seed=100 --BRT=G"
            "labrd scdz 80 70 50 100 100 100 10 --imatrix=O --seed=100 --BRT=G"
            )

set(VERIFICATION_TEST_CASES 
            "labrd scdz 75 110 50 150 150 150 10 --seed=100 --BRT=V"
            "labrd scdz 10 11 5 20 20 20 10 --seed=100 --BRT=V"
            "labrd scdz 80 50 40 100 100 100 10 --seed=100 --BRT=V"
            "labrd scdz 100 110 50 200 200 200 10 --seed=100 --BRT=V"
            "labrd scdz 110 100 75 150 150 150 10 --seed=100 --BRT=V"
            "labrd scdz 10 8 5 20 20 20 10 --seed=100 --BRT=V"
            "labrd scdz 100 180 80 200 200 200 10 --seed=100 --BRT=V"
            "labrd scdz 190 170 150 200 200 200 10 --seed=100 --BRT=V"
            "labrd scdz 200 190 150 200 200 200 10 --seed=100 --BRT=V"
            "labrd scdz 80 70 50 100 100 100 10 --seed=100 --BRT=V"
            "labrd scdz 190 170 150 200 200 200 10 --imatrix=U --seed=100 --BRT=V"
            "labrd scdz 200 190 150 200 200 200 10 --imatrix=O --seed=100 --BRT=V"
            "labrd scdz 80 70 50 100 100 100 10 --imatrix=O --seed=100 --BRT=V"
            )

set(TEST_NUM 1)
foreach(gt_test_cases IN LISTS GROUND_TRUTH_TEST_CASES)
    set(COMMANDLINE_PARAMS ${gt_test_cases})
    string(REGEX MATCH "^([^ ]+)" TEST_NAME ${gt_test_cases})
    set(TEST_NAME GROUND_TRUTH_TEST_CASE_${TEST_NUM}_${TEST_NAME})
    string(REPLACE "" ";" PROJECT_NAME_ ${PROJECT_NAME})
    add_test(${TEST_NAME} sh -c "if [ \"\$RUN_BRT_GT\" = \"TRUE\" ] || [ \"\$REPEATABILITY_TEST\" = \"TRUE\" ]; then ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME_} ${COMMANDLINE_PARAMS}; else exit 127; fi" VERBATIM)
    set_tests_properties(${TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL;No test was run, give valid arguments")
    set_tests_properties(${TEST_NAME} PROPERTIES SKIP_RETURN_CODE 127)
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "brt_gt_tests;repeatability_tests")
MATH(EXPR TEST_NUM "${TEST_NUM}+1")
endforeach()

foreach(v_test_cases IN LISTS VERIFICATION_TEST_CASES)
    set(COMMANDLINE_PARAMS ${v_test_cases})
    string(REGEX MATCH "^([^ ]+)" TEST_NAME ${v_test_cases})
    set(TEST_NAME VERIFICATION_TEST_CASE_${TEST_NUM}_${TEST_NAME})
    string(REPLACE "" ";" PROJECT_NAME_ ${PROJECT_NAME})
    add_test(${TEST_NAME} sh -c "if [ \"\$RUN_BRT_V\" = \"TRUE\" ] || [ \"\$REPEATABILITY_TEST\" = \"TRUE\" ]; then ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME_} ${COMMANDLINE_PARAMS}; else exit 127; fi" VERBATIM)
    set_tests_properties(${TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL;No test was run, give valid arguments")
    set_tests_properties(${TEST_NAME} PROPERTIES SKIP_RETURN_CODE 127)
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "brt_v_tests;repeatability_tests")
MATH(EXPR TEST_NUM "${TEST_NUM}+1")
endforeach()