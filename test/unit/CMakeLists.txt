add_executable(unitTest
        autocorr_test.cpp
        clk_test.cpp
        mmc_config_test.cpp
        newsud_test.cpp
        pseudorandom_test.cpp
        recombo_criteria_test.cpp
        writhe_test.cpp)
target_link_libraries(unitTest recombo_core gtest gtest_main)
target_compile_definitions(unitTest PRIVATE
        TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data")
target_include_directories(unitTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
include(GoogleTest)
gtest_discover_tests(unitTest)
