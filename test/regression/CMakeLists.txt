add_executable(regressionTest
        bfacf3_test.cpp
        bfacf_critical_z_test.cpp
        bfacf_regression_test.cpp
        bfacf_step_test.cpp
        bfacf_stepq_test.cpp
        bfacf_topology_invariance_test.cpp)
target_link_libraries(regressionTest recombo_core gtest gtest_main)
target_compile_definitions(regressionTest PRIVATE
        TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data"
        REPO_ROOT_DIR="${CMAKE_SOURCE_DIR}")
target_include_directories(regressionTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
include(GoogleTest)
gtest_discover_tests(regressionTest)
