add_executable(regressionTest
        bfacf_regression_test.cpp)
target_link_libraries(regressionTest recombo_core gtest gtest_main)
include(GoogleTest)
gtest_discover_tests(regressionTest)
