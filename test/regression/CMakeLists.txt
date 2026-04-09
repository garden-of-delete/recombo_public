add_executable(regressionTest
        bfacf_regression_test.cpp)
target_link_libraries(regressionTest recombo_core gtest gtest_main)
add_test(NAME regressionTest COMMAND $<TARGET_FILE:regressionTest>)
