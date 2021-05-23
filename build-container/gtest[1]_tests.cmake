add_test( Step3Tester.MakeGrid /workspaces/deal.ii-course-practice/build-container/gtest [==[--gtest_filter=Step3Tester.MakeGrid]==] --gtest_also_run_disabled_tests)
set_tests_properties( Step3Tester.MakeGrid PROPERTIES WORKING_DIRECTORY /workspaces/deal.ii-course-practice/build-container)
set( gtest_TESTS Step3Tester.MakeGrid)
